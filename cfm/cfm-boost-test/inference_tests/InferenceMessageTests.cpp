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
        //BOOST_TEST((double)msg.getIdx(5) > -A_BIG_DBL);
        BOOST_CHECK_CLOSE_FRACTION((double)msg.getIdx(12), -2.8000684380957894, tol);
        BOOST_CHECK_CLOSE_FRACTION((double)msg.getIdx(100), -1.4500684380957891, tol);
        BOOST_CHECK_CLOSE_FRACTION((double)msg.getIdx(123), -0.35006843809578925 , tol);
    }
    BOOST_AUTO_TEST_CASE(Iteration) {

        double tol = 1e-8;
        Message::const_iterator it = msg.begin();

        BOOST_CHECK_CLOSE_FRACTION((double)*it , -2.8000684380957894, tol);
        BOOST_CHECK_EQUAL(it.index(), 12);
        it ++;

        BOOST_CHECK_CLOSE_FRACTION((double)*it, -1.4500684380957891, tol);
        BOOST_CHECK_EQUAL(it.index(), 100);
        it ++;

        BOOST_CHECK_CLOSE_FRACTION((double)*it, -0.35006843809578925 , tol);
        BOOST_CHECK_EQUAL(it.index(), 123);
        it ++;

        BOOST_CHECK(it == msg.end());
    }

    BOOST_AUTO_TEST_CASE(CopyAssginment) {

        //Copy it
        Message msg2 = msg;

        //Change the old message
        msg.addToIdx(5, 0.25);

        double tol = 1e-8;

        //Check index 5
        BOOST_CHECK_CLOSE_FRACTION((double)msg.getIdx(5),  -1.910282456760932, tol);

        //Check index 123
        BOOST_CHECK_CLOSE_FRACTION((double)msg.getIdx(123), -0.510282456760932, tol);
        BOOST_CHECK_CLOSE_FRACTION((double)msg2.getIdx(123), -0.350068438095789, tol);
    }

BOOST_AUTO_TEST_SUITE_END()