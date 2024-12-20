/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentTreeNode.cpp
#
# Description: 	Contains Break and FragmentTreeNode classes, for use in
#				the fragment generation process.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FunctionalGroups.h"

const std::string FGRPS_PICKLE =
        "0 20 1e-008\n"
                "161\n"
                "fr_C=O	[CX3]=[OX1]\n"
                "fr_C=O_noCOO	[C!$(C-[OH])]=O\n"
                "fr_Al_OH	[C!$(C=O)]-[OH]\n"
                "fr_Ar_OH	c[OH1]\n"
                "fr_methoxy	[OX2](-[#6])-[CH3]\n"
                "fr_oxime	[CX3]=[NX2]-[OX2]\n"
                "fr_ester	[#6][CX3](=O)[OX2H0][#6]\n"
                "fr_Al_COO	C-C(=O)[O;H1,-]\n"
                "fr_Ar_COO	c-C(=O)[O;H1,-]\n"
                "fr_COO	[#6]C(=O)[O;H,-1]\n"
                "fr_COO2	[CX3](=O)[OX1H0-,OX2H1]\n"
                "fr_ketone	[#6][CX3](=O)[#6]\n"
                "fr_ether	[OD2]([#6])[#6]\n"
                "fr_phenol	[OX2H]-c1ccccc1\n"
                "fr_aldehyde	[CX3H1](=O)[#6]\n"
                "fr_quatN	[$([NX4+]),$([NX4]=*)]\n"
                "fr_NH2	[NH2,nH2]\n"
                "fr_NH1	[NH1,nH1]\n"
                "fr_NH0	[NH0,nH0]\n"
                "fr_Ar_N	n\n"
                "fr_Ar_NH	[nH]\n"
                "fr_aniline	c-[NX3]\n"
                "fr_Imine	[Nv3](=C)-[#6]\n"
                "fr_nitrile	[NX1]#[CX2]\n"
                "fr_hdrzine	[NX3]-[NX3]\n"
                "fr_hdrzone	C=N-[NX3]\n"
                "fr_nitroso	[N!$(N-O)]=O\n"
                "fr_N-O	[N!$(N=O)](-O)-C\n"
                "fr_nitro	[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]\n"
                "fr_azo	[#6]-N=N-[#6]\n"
                "fr_diazo	[N+]#N\n"
                "fr_azide	[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]\n"
                "fr_amide	C(=O)-N\n"
                "fr_priamide	C(=O)-[NH2]\n"
                "fr_amidine	C(=N)(-N)-[!#7]\n"
                "fr_guanido	C(=N)(N)N\n"
                "fr_imide	N(-C(=O))-C=O\n"
                "fr_isocyan	N=C=O\n"
                "fr_isothiocyan	N=C=S\n"
                "fr_thiocyan	S-C#N\n"
                "fr_halogen	[#9,#17,#35,#53]\n"
                "fr_alkyl_halide	[CX4]-[Cl,Br,I,F]\n"
                "fr_sulfide	[SX2](-[#6])-C\n"
                "fr_SH	[SH]\n"
                "fr_C=S	C=[SX1]\n"
                "fr_sulfone	S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]\n"
                "fr_sulfonamd	N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]\n"
                "fr_barbitur	C1C(=O)NC(=O)NC1=O\n"
                "fr_urea	C(=O)(-N)-N\n"
                "fr_term_acetylene	C#[CH]\n"
                "fr_imidazole	n1cncc1\n"
                "fr_furan	o1cccc1\n"
                "fr_thiophene	s1cccc1\n"
                "fr_thiazole	c1scnc1\n"
                "fr_oxazole	c1ocnc1\n"
                "fr_pyridine	n1ccccc1\n"
                "fr_piperdine	N1CCCCC1\n"
                "fr_piperzine	N1CCNCC1\n"
                "fr_morpholine	O1CCNCC1\n"
                "fr_lactam	N1C(=O)CC1\n"
                "fr_lactone	[C&R1](=O)[O&R1][C&R1]\n"
                "fr_tetrazole	c1nnnn1\n"
                "fr_epoxide	O1CC1\n"
                "fr_unbrch_alkane	[R0;D2][R0;D2][R0;D2][R0;D2]\n"
                "fr_bicyclic	[R2][R2]\n"
                "fr_benzene	c1ccccc1\n"
                "fr_phos_acid	[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]\n"
                "fr_phos_ester	[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]\n"
                "fr_nitro_arom	[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1)]\n"
                "fr_nitro_arom_nonortho	[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1);!$(cc-!:*)]\n"
                "fr_dihydropyridine	[$([NX3H1]1-C=C-C-C=C1),$([Nv3]1=C-C-C=C-C1),$([Nv3]1=C-C=C-C-C1),$([NX3H1]1-C-C=C-C=C1)]\n"
                "fr_phenol_noOrthoHbond	[$(c1(-[OX2H])ccccc1);!$(cc-!:[CH2]-[OX2H]);!$(cc-!:C(=O)[O;H1,-]);!$(cc-!:C(=O)-[NH2])]\n"
                "fr_Al_OH_noTert	[$(C-[OX2H]);!$([CX3](-[OX2H])=[OX1]);!$([CD4]-[OX2H])]\n"
                "fr_benzodiazepine	[c&R2]12[c&R1][c&R1][c&R1][c&R1][c&R2]1[N&R1][C&R1][C&R1][N&R1]=[C&R1]2\n"
                "fr_para_hydroxylation	[$([cH]1[cH]cc(c[cH]1)~[$([#8,$([#8]~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)~[$([#7X3,$([#7](~[H,c,C])~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)-!:[$([NX3H,$(NC(=O)[H,c,C])])])]\n"
                "fr_allylic_oxid	[$(C=C-C);!$(C=C-C-[N,O,S]);!$(C=C-C-C-[N,O]);!$(C12=CC(=O)CCC1C3C(C4C(CCC4)CC3)CC2)]\n"
                "fr_aryl_methyl	[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]\n"
                "fr_Ndealkylation1	[$(N(-[CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)]),$(N(-[CH2][CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)])]\n"
                "fr_Ndealkylation2	[$([N&R1]1(-C)CCC1),$([N&R1]1(-C)CCCC1),$([N&R1]1(-C)CCCCC1),$([N&R1]1(-C)CCCCCC1),$([N&R1]1(-C)CCCCCCC1)]\n"
                "fr_alkyl_carbamate	C[NH1]C(=O)OC\n"
                "fr_ketone_Topliss	[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]\n"
                "fr_ArN	[$(a-[NX3H2]),$(a-[NH1][NH2]),$(a-C(=[OX1])[NH1][NH2]),$(a-C(=[NH])[NH2])]\n"
                "fr_HOCCN	[$([OX2H1][CX4][CX4H2][NX3&R1]),$([OH1][CX4][CX4H2][NX3][CX4](C)(C)C)]\n"
                "y_2	[#6][OX2H1]\n"
                "y_3	[*;a]\n"
                "y_4	[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]\n"
                "y_5	[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]\n"
                "y_7	[#6;X3H1](-[#6;a])=[O;v2X1]\n"
                "y_9	[c;R1]1[c;R1][n;R1][c;R1][n;R1]1\n"
                "y_11	[Cl;X1][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]\n"
                "y_12	C#N\n"
                "y_13	[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]\n"
                "y_14	[#8;A;X1-,X2H1,X2C][P;X4]([#8;A;X1-,X2H1,X2C])(=[O;X1])c:c\n"
                "y_15	[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]\n"
                "y_16	[CX2]#[CX2]\n"
                "y_17	[CX3]=[CX2]=[CX3]\n"
                "y_19	[OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]\n"
                "y_20	[SX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]\n"
                "y_21	[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]\n"
                "y_22	[c][OX2][c]\n"
                "y_23	[SX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]\n"
                "y_24	[c][SX2][c]\n"
                "y_25	[#6;A;X4H1][#8;A;H1X2]\n"
                "y_28	[CX3;$([R0][#6]),$([H1R0])](=[OX1])[NX3H2]\n"
                "y_29	[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]\n"
                "y_30	[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]\n"
                "y_31	[F,Cl,Br,I][c]\n"
                "y_32	[#6]~[#7,#8,#16]\n"
                "y_33	[*,#1;*,#1;O,C,#1,N,X]-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3](-[*,#1;*,#1;O,C,#1,N,X])=[O;X1]\n"
                "y_36	[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]-1=[O;X1]\n"
                "y_37	[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1](=O)-[#6;R1]=[#6;R1]-1\n"
                "y_38	[#6;X4][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]\n"
                "y_39	[#6;a][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]\n"
                "y_42	[#6]S[#6]\n"
                "y_43	[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]\n"
                "y_44	[$([SX3](=[OX1])([#6])[#6]),$([SX3+]([OX1-])([#6])[#6])]\n"
                "y_45	[#6][S;X3]([#16;X2])=[O;X1]\n"
                "y_46	[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]\n"
                "y_48	[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]\n"
                "y_50	[#6;a]-[#7H1]-[#7H1]-[#6;a]\n"
                "y_51	[#6]-[#16;H1X2]\n"
                "y_52	[SX2H][CX4;!$(C([SX2H])~[O,S,#7,#15])]\n"
                "y_53	[C;X4H2]([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]\n"
                "y_54	[C;X4H1]([#6;H3X4])([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]\n"
                "y_55	[#6]-[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]\n"
                "y_57	[OX2r3]1[#6r3][#6r3]1\n"
                "y_58	[#6;a]-[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]\n"
                "y_59	[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1\n"
                "y_60	[c;R1]1[c;R1][n;R1][c;R1][n;R1][c;R1]1\n"
                "y_62	[SX2][CX2]#[NX1]\n"
                "y_63	[CX3]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]\n"
                "y_64	[CX3]=[CX2]=[OX1]\n"
                "y_65	[#8;H0X2]-[#6]([#6,#1;A])-[#6]=[#6]\n"
                "y_72	[NX3H2+0,NX4H3+]c\n"
                "y_75	[OX2D2][OX2D2]\n"
                "y_76	[NX2](=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])])[OX2][#6;!$(C=[#7,#8])]\n"
                "y_78	[#8;X1-][S;X3](=[O;X1])[$([CX3;!R][#6]),$([CX3H;!R])]=[#7;NX2H1,$([#7][#6;!$(C=[O,N,S])])]\n"
                "y_79	[#6;a][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1]\n"
                "y_80	[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]\n"
                "y_81	[OX2H][CX4H1;!$(C([OX2H])[O,S,#7,#15])]\n"
                "y_82	[OX2H][CX4D4;!$(C([OX2H])[O,S,#7,#15])]\n"
                "y_84	[#6]-[#7;X2]=[O;X1]\n"
                "y_85	[#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]\n"
                "y_86	[#16;X2H1]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]\n"
                "y_87	[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4;$([H1]),$([H0][#6])](=[O;X1])=[O;X1]\n"
                "y_88	[#6;R1]-,=1-,=[#6;R1]-,=[#6]-,=2-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=2-,=[#6;R1]-,=1\n"
                "y_89	N1CCNC1\n"
                "y_90	[$([SX3](=[OX1])([#6])[#6;a]),$([SX3+]([OX1-])([#6])[#6;a])]\n"
                "y_91	[$([SX4](=[OX1])(=[OX1])([#6])[#6;a])]\n"
                "y_92	C1=CSC=C1\n"
                "y_94	S1C=CN=C1\n"
                "y_98	C1OCOC1\n"
                "y_102	[#6]~[#16]\n"
                "y_103	[#7;v3X3H1]([#6;A;X4])-[$([cX3](:*):*),$([cX2+](:*):*)]\n"
                "y_104	[$([#6,#7,#8;A;-][N+]#C),$([#6-]=[N+]=[#6,#7,#8;A]),$([#6,#7,#8;A;+][#7]=[#6-]),$([#6]-[#7]=[#6,#7,#8;A])].[#6]\n"
                "y_107	[a;!c]\n"
                "-X	*-[#9,#17,#35,#53]\n"
                "-tBu	*-[C;D4]([C;D1])([C;D1])-[C;D1]\n"
                "-CF3	*-[C;D4](F)(F)F\n"
                "-C#CH	*-[C;D2]#[C;D1;H]\n"
                "-cPropyl	*-[C;D3]1-[C;D2]-[C;D2]1\n";

const std::string EXTRA_FGRPS_PICKLE =
        "0 20 1e-008\n"
                "13\n"
                "x_CH0	[CH0]\n"
                "x_CH1	[CH1]\n"
                "x_CH2	[CH2]\n"
                "x_CH3	[CH3]\n"
                "x_cH0	[cH0]\n"
                "x_cH1	[cH1]\n"
                "x_r3	[r3]\n"
                "x_r4	[r4]\n"
                "x_r5	[r5]\n"
                "x_r6	[r6]\n"
                "x_r7	[r7]\n"
                "x_r8	[r8]\n"
                "x_N4	[N+H0]\n";

const std::string PI_BOND_FGRPS_PICKLE =
        "0 20 1e-008\n"
        "1\n"
        "NitroGroup	[N+](=O)[O-]\n";
        //"azide_ion	[$([#7-]=[N+]=[#7-])]\n";
