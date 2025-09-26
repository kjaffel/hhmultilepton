# coding: utf-8

"""
HH -> multi-leptons selection methods.
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import

ak = maybe_import("awkward")


#
# dummy selector
#

@categorizer(uses={"event"})
def cat_all(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # keep all events
    return events, ak.ones_like(events.event) == 1

#
# lepton channels
#

@categorizer(uses={"channel_id"})
def cat_etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.etau.id


@categorizer(uses={"channel_id"})
def cat_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.mutau.id


@categorizer(uses={"channel_id"})
def cat_tautau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.tautau.id


@categorizer(uses={"channel_id"})
def cat_ee(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.ee.id


@categorizer(uses={"channel_id"})
def cat_mumu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.mumu.id


@categorizer(uses={"channel_id"})
def cat_emu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.emu.id

# multilepton


@categorizer(uses={"channel_id"})
def cat_c3e(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.c3e.id


@categorizer(uses={"channel_id"})
def cat_c2emu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.c2emu.id


@categorizer(uses={"channel_id"})
def cat_ce2mu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.ce2mu.id


@categorizer(uses={"channel_id"})
def cat_c3mu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.c3mu.id


@categorizer(uses={"channel_id"})
def cat_c4e(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.c4e.id


@categorizer(uses={"channel_id"})
def cat_c3emu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.c3emu.id


@categorizer(uses={"channel_id"})
def cat_c2e2mu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.c2e2mu.id


@categorizer(uses={"channel_id"})
def cat_ce3mu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.ce3mu.id


@categorizer(uses={"channel_id"})
def cat_c4mu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.ce3mu.id


# to be implemented
# @categorizer(uses={"channel_id"})
# def cat_c3etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.c3etau.id

# @categorizer(uses={"channel_id"})
# def cat_c2emutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.c2emutau.id

# @categorizer(uses={"channel_id"})
# def cat_ce2mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.ce2mutau.id

# @categorizer(uses={"channel_id"})
# def cat_c3mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.c3mutau.id

# @categorizer(uses={"channel_id"})
# def cat_c2e2tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.c2e2tau.id

# @categorizer(uses={"channel_id"})
# def cat_cemu2tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.cemu2tau.id

# @categorizer(uses={"channel_id"})
# def cat_c2mu2tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.c2mu2tau.id

# @categorizer(uses={"channel_id"})
# def cat_ce3tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.ce3tau.id

# @categorizer(uses={"channel_id"})
# def cat_cmu3tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.cmu3tau.id

# @categorizer(uses={"channel_id"})
# def cat_c4tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#     return events, events.channel_id == self.config_inst.channels.n.c4tau.id


# 3l/4l inclusive, later split into CR / SR via Z-peak
@categorizer(uses={"channel_id"})
def cat_3l0tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    catmask = events.channel_id == self.config_inst.channels.n.c3e.id
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c3mu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c2emu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.ce2mu.id)
    return events, catmask


@categorizer(uses={"channel_id"})
def cat_4l(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    catmask = events.channel_id == self.config_inst.channels.n.c4e.id
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c3emu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c2e2mu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.ce3mu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c4mu.id)
    return events, catmask

#
# QCD regions
#

@categorizer(uses={"leptons_os"})
def cat_os(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # oppositive sign leptons
    return events, events.leptons_os == 1


@categorizer(uses={"leptons_os"})
def cat_ss(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # same sign leptons
    return events, events.leptons_os == 0


@categorizer(uses={"tau2_isolated"})
def cat_iso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # isolated tau2
    return events, events.tau2_isolated == 1


@categorizer(uses={"tau2_isolated"})
def cat_noniso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # noon-isolated tau2
    return events, events.tau2_isolated == 0


#
# kinematic regions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1


@categorizer(uses={"Jet.{pt,phi}"})
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # two or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 2


@categorizer(uses={"Jet.btagPNetB"})
def cat_res1b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # exactly pnet b-tags
    wp = self.config_inst.x.btag_working_points["particleNet"]["medium"]
    tagged = events.Jet.btagPNetB > wp
    return events, ak.sum(tagged, axis=1) == 1


@categorizer(uses={"Jet.btagPNetB"})
def cat_res2b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # at least two medium pnet b-tags
    wp = self.config_inst.x.btag_working_points["particleNet"]["medium"]
    tagged = events.Jet.btagPNetB > wp
    return events, ak.sum(tagged, axis=1) >= 2


@categorizer(uses={cat_res1b, cat_res2b, "FatJet.{pt,phi}"})
def cat_boosted(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # not res1b or res2b, and exactly one selected fat jet that should also pass a tighter pt cut
    # note: this is just a draft
    mask = (
        (ak.num(events.FatJet, axis=1) == 1) &
        (ak.sum(events.FatJet.pt > 350, axis=1) == 1) &
        ~self[cat_res1b](events, **kwargs)[1] &
        ~self[cat_res2b](events, **kwargs)[1]
    )
    return events, mask


@categorizer(uses={"{Electron,Muon,Tau}.{pt,eta,phi,mass}"})
def cat_dy(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # e/mu driven DY region: mll > 40 and met < 30 (to supress tau decays into e/mu)
    leps = ak.concatenate([events.Electron * 1, events.Muon * 1, events.Tau * 1], axis=1)[:, :2]
    mask = (
        (leps.sum(axis=1).mass > 40) &
        (events[self.config_inst.x.met_name].pt < 30)
    )
    return events, mask


@cat_dy.init
def cat_dy_init(self: Categorizer) -> None:
    self.uses.add(f"{self.config_inst.x.met_name}.{{pt,phi}}")


@categorizer(uses={"{Electron,Muon,Tau}.{pt,eta,phi,mass}"})
def cat_tt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # tt region: met > 30 (due to neutrino presence in leptonic w decays)
    mask = events[self.config_inst.x.met_name].pt > 30
    return events, mask


@cat_tt.init
def cat_tt_init(self: Categorizer) -> None:
    self.uses.add(f"{self.config_inst.x.met_name}.{{pt,phi}}")
