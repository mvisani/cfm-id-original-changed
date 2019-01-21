/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for features.cpp
#
# Author: Felicity Allen, Fei Wang
# Created: November 2012, 2019
#########################################################################*/
#include <boost/test/unit_test.hpp>

#include "Message.h"

struct MessageFixture{
    MessageFixture() {
        msg = Message(150);
        msg.addToIdx(100, 0.55);
        msg.addToIdx(123, 1.65);
        msg.addToIdx(12, -0.8);
    }
    ~MessageFixture() {};
    Message msg;
};

BOOST_FIXTURE_TEST_SUITE(MessageTests, MessageFixture)

    BOOST_AUTO_TEST_CASE(Basics) {
        double tol = 1e-8;
        // TO DO FIX
        //BOOST_CHECK_CLOSE(std::exp(msg.getIdx(5)), -A_BIG_DBL , tol);
        BOOST_CHECK_CLOSE((double)msg.getIdx(12), -2.8000684380957894, tol);
        BOOST_CHECK_CLOSE((double)msg.getIdx(100), -1.4500684380957891, tol);
        BOOST_CHECK_CLOSE((double)msg.getIdx(123), -0.35006843809578925 , tol);
    }
    BOOST_AUTO_TEST_CASE(Iteration) {

        double tol = 1e-8;
        Message::const_iterator it = msg.begin();

        BOOST_CHECK_CLOSE((double)*it , -2.8000684380957894, tol);
        BOOST_CHECK_EQUAL(it.index(), 12);
        it ++;

        BOOST_CHECK_CLOSE((double)*it, -1.4500684380957891, tol);
        BOOST_CHECK_EQUAL(it.index(), 100);
        it ++;

        BOOST_CHECK_CLOSE((double)*it, -0.35006843809578925 , tol);
        BOOST_CHECK_EQUAL(it.index(), 123);
        it ++;

        BOOST_CHECK(it == msg.end());
    }

    BOOST_AUTO_TEST_CASE(CopyAssginment) {

        //Copy it
        Message msg2 = msg;

        //Change the old message
        msg.addToIdx(5, 0.25);

        //Check index 5
        /*double tmp = m.getIdx(5);
        double tmp2 = m2.getIdgiux(5);
        if( fabs( tmp - tmp2 ) < tol || tmp2 > -A_BIG_DBL || fabs( tmp - -1.910282456760932) > tol){
            std::cout << "Unexpected values at index 5 after copy: " << tmp << " vs " << tmp2 << std::endl;
            pass = false;
        }

        //Check index 123
        tmp = m.getIdx(123);
        tmp2 = m2.getIdx(123);
        if( fabs( tmp - -0.510282456760932) > tol || fabs( tmp2 - -0.350068438095789) > tol ){
            std::cout << "Unexpected values at index 123 after copy: " << tmp << " vs " << tmp2 << std::endl;
            pass = false;
        }
        */
    }

BOOST_AUTO_TEST_SUITE_END()