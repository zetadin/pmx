#!/usr/bin/env python
from pmx.estimators import BAR, Jarz, Crooks
import pickle
from numpy.testing import assert_almost_equal
from numpy.random import seed


def test_BAR(gf):
    # these pickled files had to be created with py3 using protocol=2
    # in order to be py2/3 compatible
    wf = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))
    wr = pickle.load(open(gf("dgdl/wr.pkl"), "rb"))

    seed(42)
    bar = BAR(wf=wf, wr=wr, T=298, nblocks=3, nboots=10)
    assert_almost_equal(bar.dg, 0.8532451987153411, decimal=7)
    assert_almost_equal(bar.err, 0.557805610661955, decimal=7)
    assert_almost_equal(bar.conv, -0.1184253367468773, decimal=7)
    assert_almost_equal(bar.err_blocks, 1.2327258529433998, decimal=7)
    assert_almost_equal(bar.err_boot, 0.4613033210812272, decimal=7)
    assert_almost_equal(bar.conv_err_boot, 0.05189927533062071, decimal=7)


def test_Jarz(gf):
    wf = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))
    wr = pickle.load(open(gf("dgdl/wr.pkl"), "rb"))

    seed(42)
    jarz = Jarz(wf=wf, wr=wr, T=298, nblocks=3, nboots=10)
    assert_almost_equal(jarz.dg_mean, 0.47552516706599357, decimal=7)
    assert_almost_equal(jarz.dg_for, -0.036297126950559477, decimal=7)
    assert_almost_equal(jarz.dg_rev, 0.98734746108254667, decimal=7)
    assert_almost_equal(jarz.err_blocks_for, 1.5034007401843161, decimal=7)
    assert_almost_equal(jarz.err_blocks_rev, 0.8732303661985571, decimal=7)
    assert_almost_equal(jarz.err_boot_for, 1.31509469126237, decimal=7)
    assert_almost_equal(jarz.err_boot_rev, 0.6698096954291023, decimal=7)


def test_Crooks(gf):
    wf = pickle.load(open(gf("dgdl/wf.pkl"), "rb"))
    wr = pickle.load(open(gf("dgdl/wr.pkl"), "rb"))

    seed(42)
    cgi = Crooks(wf=wf, wr=wr, nblocks=3, nboots=10)
    assert_almost_equal(cgi.dg, 0.9312390509679932, decimal=7)
    assert_almost_equal(cgi.err_blocks, 1.4890168005551732, decimal=7)
    assert_almost_equal(cgi.err_boot1, 0.6492475842858807, decimal=7)
    assert_almost_equal(cgi.err_boot2, 0.8283881591989153, decimal=7)
