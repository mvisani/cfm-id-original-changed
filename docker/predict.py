import argparse
import concurrent.futures
import logging
import multiprocessing
import os
import subprocess
import time
import re

from collections import deque
from typing import Tuple, Dict

class AdductType(object):
    M_H_POS = '[M+H]+'
    M_H_NEG = '[M-H]-'
    ALL = 'ALL'

class IonizationType(object):
    EI = 'EI'
    ESI = 'ESI'

class ElapsedFormatter(logging.Formatter):
    def __init__(self):
        super().__init__()
        self.start_time = time.time()

    def format(self, record):
        # using timedelta here for convenient default formatting
        time_str = self.get_time_between(record.created)
        return "[{}] [{}] {}".format(time_str, record.levelname, record.getMessage())

    def get_time_between(self, current):
        hours, rem = divmod(current - self.start_time, 3600)
        minutes, seconds = divmod(rem, 60)
        return "{:0>2}:{:0>2}:{:0>2}".format(int(hours), int(minutes), int(seconds))


def esi_prediction_task(smiles, output_id, task_config) -> Tuple[
    float, str, Dict]:

    # record start time
    msrb_predicted_adduct_types = []
    msml_predicted_adduct_types = []
    msml_adduct_types = []

    msrb_binary = task_config['msrb_binary']
    need_msrb = task_config['need_msrb']
    need_msml = task_config['need_msml']
    adduct_type = task_config['adduct_type']
    output_to_file = task_config['output_to_file']

    predicted_spectra_str = ''
    if need_msrb:
        # if msrb only mode or msml + msrb
        # Create msrb command
        msrb_cmd = ['java', '-jar', msrb_binary,'-ismi', smiles]

        if output_to_file:
            msrb_cmd.append('-presult')
            msrb_cmd.append(output_id)
        else:
            msrb_cmd.append('-o')
            msrb_cmd.append(output_id)

        if adduct_type != AdductType.ALL:
            msrb_cmd.append('-a')
            msrb_cmd.append(adduct_type)
            msrb_cmd.append('-na')

        # run msrb first
        # print(' '.join(msrb_cmd))
        process = subprocess.Popen(
            msrb_cmd, stdout=subprocess.PIPE, shell = False)
        std_output, std_error = process.communicate()

        if process.returncode != 0:
            print('Error with MSRB: ', ' '.join(msrb_cmd))
            print(std_error.decode('ascii'))
        
        msrb_results_rows = [row for row in std_output.decode('ascii').split('\n') if '[main]' not in row and 'STATUS REPORT' not in row and 'chemical class' not in row] 
        msrb_results_str = '\n'.join(msrb_results_rows)
        predicted_spectra_str = msrb_results_str
        
        if  adduct_type in [AdductType.ALL, AdductType.M_H_POS]:
            if AdductType.M_H_NEG not in msrb_results_str:
                msml_adduct_types.append(AdductType.M_H_POS)
            else:
                msrb_predicted_adduct_types.append(AdductType.M_H_POS)

        if  adduct_type in [AdductType.ALL, AdductType.M_H_NEG]:
            if AdductType.M_H_NEG not in msrb_results_str:
                msml_adduct_types.append(AdductType.M_H_NEG)
            else:
                msrb_predicted_adduct_types.append(AdductType.M_H_POS)

    elif adduct_type == AdductType.ALL: 
        # we don't need msrb and we need all adduct msml can do
        msml_adduct_types = [AdductType.M_H_POS, AdductType.M_H_NEG]
    else:
        msml_adduct_types = [adduct_type]

    # Run msml
    if need_msml:
        include_annotations =  task_config['cfm_include_annotations']
        postprocessing_method = str(task_config['cfm_postprocessing_method'])
        suppress_exception = '1'
        postprocessing_energy = str(task_config['cfm_postprocessing_energy'])
        min_peak_intensity=  str(task_config['cfm_min_peak_intensity'])
        override_min_peaks=  str(task_config['cfm_override_min_peaks'])
        override_max_peaks=  str(task_config['cfm_override_max_peaks'])

        for adduct_type in msml_adduct_types:
            output_file = '{}_{}.log'.format(output_file, adduct_type) if output_to_file is False else 'stdout'
            msml_cmd = [task_config['cfm_predict_binary'], smiles, '0.001', task_config[IonizationType.ESI][adduct_type]['param_file'],
                        task_config[IonizationType.ESI][adduct_type]['config_file'],include_annotations, output_file,
                        postprocessing_method, suppress_exception, postprocessing_energy, 
                        min_peak_intensity, override_min_peaks, override_max_peaks, output_id]
            process = subprocess.Popen(
                msml_cmd, stdout=subprocess.PIPE, shell=os.name == 'nt', cwd='./')
            
            #TODO Catch Error and send as output
            std_output, std_error = process.communicate()

            if process.returncode != 0:
                print('Error with MSML: ', ' '.join(msml_cmd))
                print(std_error.decode('ascii'))

            if not output_to_file and std_error == None:
                predicted_spectra_str += std_output.decode('ascii')
                msml_predicted_adduct_types.append(adduct_type)
    
    # replace null id for msrb and msml
    if not output_to_file:
        predicted_spectra_str = predicted_spectra_str.replace("#ID=null", "#ID={}".format(output_id))
        predicted_spectra_str = predicted_spectra_str.replace("#ID=NullId", "#ID={}".format(output_id))

    return predicted_spectra_str
    
if __name__ == '__main__':

    # get arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input_list', help='list of input compounds', default='./public/input.txt')
    parser.add_argument('--adduct_type', help='one of [ALL|[M+H]+|[M-H]-]', type=str, default='ALL')
    parser.add_argument('--output', help='one of [output dir|output file(ends with .txt or .log)| stdout]', type=str, default='stdout')
    
    # CFM MSML configs
    parser.add_argument('--cfm_predict_binary',
                        help='location of cfm-predict', default='cfm-predict')
    parser.add_argument('--cfm_esi_pos_param_file',
                        help='location of cfm msml [M+H]+ model param file', default='/trained_models_cfmid4.0/ESI/[M+H]+/param_output.log')
    parser.add_argument('--cfm_esi_pos_config_file',
                        help='location of cfm msml [M+H]+ model config file', default='/trained_models_cfmid4.0/ESI/[M+H]+/param_config.txt')
    parser.add_argument('--cfm_esi_neg_param_file',
                        help='location of cfm msml [M-H]- model param file', default='/trained_models_cfmid4.0/ESI/[M-H]-/param_output.log')
    parser.add_argument('--cfm_esi_neg_config_file',
                        help='location of cfm msml [M-H]- model config file', default='/trained_models_cfmid4.0/ESI/[M-H]-/param_config.txt')

    parser.add_argument('--cfm_include_annotations', type=str, default='0')
    parser.add_argument('--cfm_postprocessing_method', help='postprocessing method', type=str, default='2')
    parser.add_argument('--cfm_postprocessing_energy', help='postprocessing energy', type=int, default=80)
    parser.add_argument('--cfm_min_peak_intensity', help='min_peak_intensity [0,100.0] ', type=int, default=0)
    parser.add_argument('--cfm_override_min_peaks', help='override_min_peaks', type=int, default=5)
    parser.add_argument('--cfm_override_max_peaks', help='override_max_peaks', type=int, default=30)

    # MSRB config
    parser.add_argument('--msrb_binary', help='location of msrb-fragmenter jar file',
                        default='./msrb-fragmenter.jar')

    args = parser.parse_args()
    # get time
    start_time = time.time()

    # get input
    input_list = args.input_list
    adduct_type = args.adduct_type

    # start or numbewr of process
    num_processes = multiprocessing.cpu_count()

    # build output folder
    is_std_output = False
    is_single_output = False
    output = os.path.abspath(args.output)
    if output == 'stdout':
        is_std_output = True
    elif '.txt' in output or '.log' in output:
        is_single_output = True
        dir_path = os.path.dirname(output)
        if dir_path != '' and not os.path.exists(dir_path):
            os.makedirs(dir_path)
    else:
        os.makedirs(output, exist_ok=True)

    log_level = logging.DEBUG
    # setup logger
    logger = logging.getLogger('cfm')
    logger.setLevel(log_level)

    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    # create formatter and add it to the handlers
    formatter = ElapsedFormatter()
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)

    # get config
    task_config = {
        'cfm_predict_binary' : args.cfm_predict_binary,
        IonizationType.ESI:{
            AdductType.M_H_POS : {'param_file' : args.cfm_esi_pos_param_file,
            'config_file' : args.cfm_esi_pos_config_file},
            AdductType.M_H_NEG : {'param_file' : args.cfm_esi_neg_param_file,
            'config_file' : args.cfm_esi_neg_config_file}
        },
        #IonizationType.EI:{
        #    'param_file' : args.cfm_ei_param_file,
        #    'config_file' : args.cfm_ei_config_file
        #},
        'cfm_include_annotations' : args.cfm_include_annotations,
        'cfm_postprocessing_method' : args.cfm_postprocessing_method,
        'cfm_postprocessing_energy' : args.cfm_postprocessing_energy,
        'cfm_min_peak_intensity': args.cfm_min_peak_intensity,
        'cfm_override_min_peaks': args.cfm_override_min_peaks, 
        'cfm_override_max_peaks': args.cfm_override_max_peaks,
        'msrb_binary': args.msrb_binary,
        'need_msml': True,
        'need_msrb': True,
        'adduct_type': adduct_type,
        'output': output,
        'output_to_file': not is_std_output and not is_single_output
    }

    # logging some info
    logger.info('=' * 120)
    logger.info('Input from: {}'.format(input_list))
    logger.info('Output to: {}'.format(output))
    logger.info('Num of Processors: {}'.format(num_processes))
    logger.info('=' * 120)
    logger.info('cfm_predict_binary: {}'.format(task_config['cfm_predict_binary']))
    logger.info('cfm_esi_pos_param_file: {}'.format(task_config[IonizationType.ESI][AdductType.M_H_POS]['param_file']))
    logger.info('cfm_esi_pos_config_file: {}'.format(task_config[IonizationType.ESI][AdductType.M_H_POS]['config_file']))
    logger.info('cfm_esi_neg_param_file: {}'.format(task_config[IonizationType.ESI][AdductType.M_H_POS]['param_file']))
    logger.info('cfm_esi_neg_config_file: {}'.format(task_config[IonizationType.ESI][AdductType.M_H_NEG]['config_file']))
    #logger.info('cfm_ei_param_file: {}'.format(task_config[IonizationType.EI]['param_file']))
    #logger.info('cfm_ei_config_file: {}'.format(task_config[IonizationType.EI]['config_file']))
    logger.info('cfm_include_annotations: {}'.format(task_config['cfm_include_annotations']))
    logger.info('cfm_postprocessing_method: {}'.format(task_config['cfm_postprocessing_method']))
    logger.info('cfm_postprocessing_energy: {}'.format(task_config['cfm_postprocessing_energy']))
    logger.info('=' * 120)
    logger.info('msrb_binary: {}'.format(args.msrb_binary))
    logger.info('=' * 120)

    finished_predict_task = 0
    finished_predict_task_last_interval = 0
    overall_process_rate = 0
    current_process_rate = 0

    # We can use a with statement to ensure threads are cleaned up promptly
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_processes) as executor, open(input_list) as file_in:
        # start a future for a thread which sends work in through the queue
        # futures_in_working = {executor.submit(get_tasks, db_file, worker_grp_id) : ""}

        futures_in_working = {}
        for row in file_in:
            items = re.split(r'[;|,|\t|\s]', row.rstrip())

            output_id,smiles = items[0], items[1]
            futures_in_working[executor.submit(
                    esi_prediction_task, smiles, output_id, task_config)] = output_id
        total_tasks = len(futures_in_working)

        while futures_in_working:

            # check for status of the futures which are currently working
            completed, not_completed = concurrent.futures.wait(
                futures_in_working, timeout=0.25,
                return_when=concurrent.futures.FIRST_COMPLETED)

            # process any completed futures
            finished_tasks = []
            predictied_spectra_record = []
            finished_predict_task_last_interval = 0
            task_used_time_total = 0

            single_output_out = open(output, 'a+') if is_single_output else None
            for future in completed:
                try:
                    data = future.result()
                except Exception as exc:
                    print('Generated an exception:', exc)
                else:
                    finished_task_id = futures_in_working[future]
                    predicted_spectra_str = data
                    finished_tasks.append(finished_task_id)
                    
                    if is_single_output:
                        single_output_out.write(predicted_spectra_str)
                        single_output_out.write('\n'*2)
                    elif is_std_output:
                        print(predicted_spectra_str)
                        print('\n'*2)

                    finished_predict_task += 1

                # remove the now completed future
                del futures_in_working[future]

            if is_single_output:
                single_output_out.close()

            # compute some metrics
            sleep_interval = 30
            if not is_std_output:
                num_not_completed = len(not_completed)
                logger.info('{} tasks in the queue. Finished {} prediction tasks since start.'.format(num_not_completed, finished_predict_task))
                logger.info('Manager sleep for {:.2f} s'.format(sleep_interval))
                time.sleep(sleep_interval)
