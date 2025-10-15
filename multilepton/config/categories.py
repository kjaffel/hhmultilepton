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
    # multilepton channels
    _add_category(config, name="3e", id=14, selection="cat_3e", label=config.channels.n.3e.label)
    _add_category(config, name="2emu", id=15, selection="cat_2emu", label=config.channels.n.2emu.label)
    _add_category(config, name="e2mu", id=16, selection="cat_e2mu", label=config.channels.n.e2mu.label)
    _add_category(config, name="3mu", id=17, selection="cat_3mu", label=config.channels.n.3mu.label)
    _add_category(config, name="4e", id=18, selection="cat_4e", label=config.channels.n.4e.label)
    _add_category(config, name="3emu", id=19, selection="cat_3emu", label=config.channels.n.3emu.label)
    _add_category(config, name="2e2mu", id=20, selection="cat_2e2mu", label=config.channels.n.2e2mu.label)
    _add_category(config, name="e3mu", id=21, selection="cat_e3mu", label=config.channels.n.e3mu.label)
    _add_category(config, name="4mu", id=22, selection="cat_4mu", label=config.channels.n.4mu.label)
    # 3l/4l inclusive, later split into CR / SR via Z-peak
    _add_category(config, name="3l0tau", id=1001, selection="cat_3l0tau", label="$3\ell 0\tau_h$")  # noqa: W605
    _add_category(config, name="4l", id=1002, selection="cat_4l", label="$4\ell$")  # noqa: W605
    # FIXME to be implemented
    # _add_category(config, name="3etau", id=16, selection="cat_3etau", label=config.channels.n.3etau.label)
    # _add_category(config, name="2emutau", id=17, selection="cat_2emutau", label=config.channels.n.2emutau.label)
    # _add_category(config, name="e2mutau", id=18, selection="cat_e2mutau", label=config.channels.n.e2mutau.label)
    # _add_category(config, name="3mutau", id=19, selection="cat_3mutau", label=config.channels.n.3mutau.label)
    # _add_category(config, name="2e2tau", id=20, selection="cat_2e2tau", label=config.channels.n.2e2tau.label)
    # _add_category(config, name="emu2tau", id=21, selection="cat_emu2tau", label=config.channels.n.emu2tau.label)
    # _add_category(config, name="2mu2tau", id=22, selection="cat_2mu2tau", label=config.channels.n.2mu2tau.label)
    # _add_category(config, name="e3tau", id=23, selection="cat_e3tau", label=config.channels.n.e3tau.label)
    # _add_category(config, name="mu3tau", id=24, selection="cat_mu3tau", label=config.channels.n.mu3tau.label)
    # _add_category(config, name="4tau", id=25, selection="cat_4tau", label=config.channels.n.4tau.label)
    # bveto
    _add_category(config, name="bveto_on", id=30001, selection="cat_bveto_on", label="bveto on")
    _add_category(config, name="bveto_off", id=30002, selection="cat_bveto_off", label="bveto off")
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

    # Loose category for BDT trainning + tight + trigmatch
    _add_category(config, 
        name="eormu", 
        id=10000, 
        selection="cat_e_or_mu", 
        label=r"e or $\mu$", 
        tags={"eormu"}
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

    #
    # build groups
    #
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
    # control categories
    control_categories = {
        # channels first
        "channel": CategoryGroup(["ee", "mumu", "emu"], is_complete=False, has_overlap=False),
        # kinematic regions in the middle (to be extended)
        "kin": CategoryGroup(["incl", "dy", "tt"], is_complete=True, has_overlap=True),
        # relative sign last
        "sign": CategoryGroup(["os"], is_complete=False, has_overlap=False),
    }
    # Creating category combinations
    sig_sideband_categories = {
        "channel": CategoryGroup(
            ["3e", "3mu", "2emu", "e2mu", "4e", "4mu", "2e2mu", "3emu", "e3mu"],
            is_complete=True,
            has_overlap=False,
        ),
        "sel": CategoryGroup(["tight", "nontight"], is_complete=False, has_overlap=False),
        "trig": CategoryGroup(["trigmatch", "nontrigmatch"], is_complete=True, has_overlap=False),
        "vetobtag": CategoryGroup(["bveto_on", "bveto_off"], is_complete=True, has_overlap=False),
        "sign": CategoryGroup(["os", "ss"], is_complete=True, has_overlap=False),
    }
    bdt_categories = {
        "loose_ch": CategoryGroup(["eormu"], is_complete=False, has_overlap=False),
        "sel": CategoryGroup(["tight_bdt", "nontight_bdt"], is_complete=False, has_overlap=False),
        "trig": CategoryGroup(["trigmatch_bdt", "nontrigmatch_bdt"], is_complete=True, has_overlap=False),
        "vetobtag": CategoryGroup(["bveto_on", "bveto_off"], is_complete=True, has_overlap=False),
    }
    
    for cnm, cdict in {
        'main': main_categories, 
        'control': control_categories,               
        'sideband': sig_sideband_categories, 
        'bdt': bdt_categories
        }.items():
        
        add_qcd_group = False
        if cnm == 'main': 
            add_qcd_group=True
        
        create_category_combinations(
            config=config,
            categories= cdict,
            name_fn=name_fn,
            kwargs_fn=functools.partial(kwargs_fn, add_qcd_group),
            skip_fn_ctrl= skip_fn_ctrl
        )

