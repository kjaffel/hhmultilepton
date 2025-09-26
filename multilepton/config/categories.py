# coding: utf-8

"""
Definition of categories.
"""

import functools

import order as od

from columnflow.config_util import add_category, create_category_combinations, CategoryGroup
from columnflow.types import Any


def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    # root category (-1 has special meaning in cutflow)
    root_cat = add_category(config, name="all", id=-1, selection="cat_all", label="")
    _add_category = functools.partial(add_category, parent=root_cat)

    # lepton channels
    _add_category(config, name="etau", id=1, selection="cat_etau", label=config.channels.n.etau.label)
    _add_category(config, name="mutau", id=2, selection="cat_mutau", label=config.channels.n.mutau.label)
    _add_category(config, name="tautau", id=3, selection="cat_tautau", label=config.channels.n.tautau.label)
    _add_category(config, name="ee", id=4, selection="cat_ee", label=config.channels.n.ee.label)
    _add_category(config, name="mumu", id=5, selection="cat_mumu", label=config.channels.n.mumu.label)
    _add_category(config, name="emu", id=6, selection="cat_emu", label=config.channels.n.emu.label)
    # new multilepton channels
    _add_category(config, name="c3e", id=14, selection="cat_c3e", label=config.channels.n.c3e.label)
    _add_category(config, name="c2emu", id=15, selection="cat_c2emu", label=config.channels.n.c2emu.label)
    _add_category(config, name="ce2mu", id=16, selection="cat_ce2mu", label=config.channels.n.ce2mu.label)
    _add_category(config, name="c3mu", id=17, selection="cat_c3mu", label=config.channels.n.c3mu.label)
    _add_category(config, name="c4e", id=18, selection="cat_c4e", label=config.channels.n.c4e.label)
    _add_category(config, name="c3emu", id=19, selection="cat_c3emu", label=config.channels.n.c3emu.label)
    _add_category(config, name="c2e2mu", id=20, selection="cat_c2e2mu", label=config.channels.n.c2e2mu.label)
    _add_category(config, name="ce3mu", id=21, selection="cat_ce3mu", label=config.channels.n.ce3mu.label)
    _add_category(config, name="c4mu", id=22, selection="cat_c4mu", label=config.channels.n.c4mu.label)
    # to be implemented
    # _add_category(config, name="c3etau", id=16, selection="cat_c3etau", label=config.channels.n.c3etau.label)
    # _add_category(config, name="c2emutau", id=17, selection="cat_c2emutau", label=config.channels.n.c2emutau.label)
    # _add_category(config, name="ce2mutau", id=18, selection="cat_ce2mutau", label=config.channels.n.ce2mutau.label)
    # _add_category(config, name="c3mutau", id=19, selection="cat_c3mutau", label=config.channels.n.c3mutau.label)
    # _add_category(config, name="c2e2tau", id=20, selection="cat_c2e2tau", label=config.channels.n.c2e2tau.label)
    # _add_category(config, name="cemu2tau", id=21, selection="cat_cemu2tau", label=config.channels.n.cemu2tau.label)
    # _add_category(config, name="c2mu2tau", id=22, selection="cat_c2mu2tau", label=config.channels.n.c2mu2tau.label)
    # _add_category(config, name="ce3tau", id=23, selection="cat_ce3tau", label=config.channels.n.ce3tau.label)
    # _add_category(config, name="cmu3tau", id=24, selection="cat_cmu3tau", label=config.channels.n.cmu3tau.label)
    # _add_category(config, name="c4tau", id=25, selection="cat_c4tau", label=config.channels.n.c4tau.label)

    # bveto
    _add_category(config, name="bveto_on", id=30001, selection="cat_bveto_on", label="bveto on")
    _add_category(config, name="bveto_off", id=30002, selection="cat_bveto_off", label="bveto off")

    # Loose category for BDT trainning + tight + trigmatch
    _add_category(config,
        name="ceormu",
        id=10000,
        selection="cat_e_or_mu",
        label=r"e or $\mu$",
        tags={"ceormu"},
    )
    # tight/nontight
    _add_category(config,
        name="tight_bdt",
        id=11000,
        selection="cat_tight_bdt",
        label="tight",
        tags={"tight_bdt"},
    )
    _add_category(config,
        name="nontight_bdt",
        id=12000,
        selection="cat_nontight_bdt",
        label="fakeable",
        tags={"nontight_bdt"},
    )
    # trigmatch
    _add_category(
        config,
        name="trigmatch_bdt",
        id=13000,
        selection="cat_trigmatch_bdt",
        label="trigger matched",
        tags={"trigmatch_bdt"},
    )
    _add_category(config,
        name="nontrigmatch_bdt",
        id=14000,
        selection="cat_nontrigmatch_bdt",
        label="trigger unmatched",
        tags={"nontrigmatch_bdt"},
    )
    # tight/nontight
    _add_category(config,
        name="tight",
        id=10001,
        selection="cat_tight",
        label="tight",
        tags={"tight"},
    )
    _add_category(config,
        name="nontight",
        id=10002,
        selection="cat_nontight",
        label="fakeable",
        tags={"nontight"},
    )
    # trigmatch
    _add_category(config,
        name="trigmatch",
        id=10003,
        selection="cat_trigmatch",
        label="trigger matched",
        tags={"trigmatch"},
    )
    _add_category(config,
        name="nontrigmatch",
        id=10004,
        selection="cat_nontrigmatch",
        label="trigger unmatched",
        tags={"nontrigmatch"},
    )
    ######################################################

    # qcd regions
    _add_category(config, name="os", id=10, selection="cat_os", label="OS", tags={"os"})
    _add_category(config, name="ss", id=11, selection="cat_ss", label="SS", tags={"ss"})
    _add_category(config, name="iso", id=12, selection="cat_iso", label=r"iso", tags={"iso"})
    _add_category(config, name="noniso", id=13, selection="cat_noniso", label=r"non-iso", tags={"noniso"})  # noqa: E501

    # kinematic categories
    _add_category(config, name="incl", id=100, selection="cat_incl", label="inclusive")
    _add_category(config, name="2j", id=110, selection="cat_2j", label="2 jets")
    _add_category(config, name="dy", id=210, selection="cat_dy", label="DY enriched")
    _add_category(config, name="tt", id=220, selection="cat_tt", label=r"$t\bar{t}$ enriched")

    _add_category(config, name="res1b", id=300, selection="cat_res1b", label="res1b")
    _add_category(config, name="res2b", id=301, selection="cat_res2b", label="res2b")
    _add_category(config, name="boosted", id=310, selection="cat_boosted", label="boosted")

    # multilepton
    # 3l/4l inclusive, later split into CR / SR via Z-peak
    _add_category(config, name="cat3l0tau", id=1001, selection="cat_3l0tau", label="$3\ell 0\tau_h$")  # noqa: W605
    _add_category(config, name="cat4l", id=1002, selection="cat_4l", label="$4\ell$")  # noqa: W605

    #
    # build groups
    #

    def name_fn(categories: dict[str, od.Category]) -> str:
        return "__".join(cat.name for cat in categories.values() if cat)

    def kwargs_fn(categories: dict[str, od.Category], add_qcd_group: bool = True) -> dict[str, Any]:
        # build auxiliary information
        aux = {}
        if add_qcd_group:
            aux["qcd_group"] = name_fn({
                name: cat for name, cat in categories.items()
                if name not in {"sign", "tau2"}
            })
        # return the desired kwargs
        return {
            # just increment the category id
            # NOTE: for this to be deterministic, the order of the categories must no change!
            "id": "+",
            # join all tags
            "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
            # auxiliary information
            "aux": aux,
            # label
            "label": ", ".join([
                cat.label or cat.name
                for cat in categories.values()
                # if cat.name != "os"  # os is the default
            ]) or None,
        }

    # main analysis categories
    main_categories = {
        # channels first
        "channel": CategoryGroup(["etau", "mutau", "tautau"], is_complete=False, has_overlap=False),
        # kinematic regions in the middle (to be extended)
        "kin": CategoryGroup(["incl", "2j", "res1b", "res2b", "boosted"], is_complete=True, has_overlap=True),
        # qcd regions last
        "sign": CategoryGroup(["os", "ss"], is_complete=True, has_overlap=False),
        "tau2": CategoryGroup(["iso", "noniso"], is_complete=True, has_overlap=False),
    }

    create_category_combinations(
        config=config,
        categories=main_categories,
        name_fn=name_fn,
        kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=True),
    )

# ##############################################################################
    # Creating category combinations
    categories_sig_sideband = {
        "channel": CategoryGroup(
            ["c3e", "c3mu", "c2emu", "ce2mu", "c4e", "c4mu", "c2e2mu", "c3emu", "ce3mu"],
            is_complete=True,
            has_overlap=False,
        ),
        "sel": CategoryGroup(["tight", "nontight"], is_complete=False, has_overlap=False),
        "trig": CategoryGroup(["trigmatch", "nontrigmatch"], is_complete=True, has_overlap=False),
        "vetobtag": CategoryGroup(["bveto_on", "bveto_off"], is_complete=True, has_overlap=False),
        "sign": CategoryGroup(["os", "ss"], is_complete=True, has_overlap=False),
    }

    create_category_combinations(
        config=config,
        categories=categories_sig_sideband,
        name_fn=name_fn,
        kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=False),
    )

    bdt_categories = {
        "loose_ch": CategoryGroup(["ceormu"], is_complete=False, has_overlap=False),
        "sel": CategoryGroup(["tight_bdt", "nontight_bdt"], is_complete=False, has_overlap=False),
        "trig": CategoryGroup(["trigmatch_bdt", "nontrigmatch_bdt"], is_complete=True, has_overlap=False),
        "vetobtag": CategoryGroup(["bveto_on", "bveto_off"], is_complete=True, has_overlap=False),
    }

    create_category_combinations(
        config=config,
        categories=bdt_categories,
        name_fn=name_fn,
        kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=False),
    )
# #############################################################################

    # control categories
    control_categories = {
        # channels first
        "channel": CategoryGroup(["ee", "mumu", "emu"], is_complete=False, has_overlap=False),
        # kinematic regions in the middle (to be extended)
        "kin": CategoryGroup(["incl", "dy", "tt"], is_complete=True, has_overlap=True),
        # relative sign last
        "sign": CategoryGroup(["os"], is_complete=False, has_overlap=False),
    }

    def skip_fn_ctrl(categories: dict[str, od.Category]) -> bool:
        if "channel" not in categories or "kin" not in categories:
            return False
        ch_cat = categories["channel"]
        kin_cat = categories["kin"]
        # skip dy in emu
        if kin_cat.name == "dy" and ch_cat.name == "emu":
            return True
        # skip tt in ee/mumu
        if kin_cat.name == "tt" and ch_cat.name in ("ee", "mumu"):
            return True
        return False

    create_category_combinations(
        config=config,
        categories=control_categories,
        name_fn=name_fn,
        kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=False),
        skip_fn=skip_fn_ctrl,
    )
