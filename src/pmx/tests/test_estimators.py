#!/usr/bin/env python
from pmx.estimators import BAR, Jarz, Crooks
import pickle
from numpy.testing import assert_almost_equal


def test_BAR(gf):
    wf = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))
    wr = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))

    bar = BAR(wf=wf, wr=wr, T=298, nblocks=3)
    assert_almost_equal(bar.dg, 3.1858683778, decimal=7)
    assert_almost_equal(bar.err, 0.4919836061, decimal=7)
    assert_almost_equal(bar.err_blocks, 1.5970471569094371, decimal=7)


def test_Jarz(gf):
    wf = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))
    wr = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))

    jarz = Jarz(wf=wf, wr=wr, T=298, nblocks=3)
    assert_almost_equal(jarz.dg_mean, 3.9864465767089428, decimal=7)
    assert_almost_equal(jarz.dg_for, -0.036297126950559477, decimal=7)
    assert_almost_equal(jarz.dg_rev, 8.0091902803684452, decimal=7)
    assert_almost_equal(jarz.err_blocks_for, 1.5034007401843161, decimal=7)
    assert_almost_equal(jarz.err_blocks_rev, 2.3918655755627714, decimal=7)


def test_Crooks(gf):
    wf = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))
    wr = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))

    cgi = Crooks(wf=wf, wr=wr, nblocks=3)
    assert_almost_equal(cgi.dg, 3.2655948136928199, decimal=7)
    assert_almost_equal(cgi.err_blocks, 1.710028287094465, decimal=7)
