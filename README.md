# HH → Multilepton

## Introduction

This is the code base for the Run2+Run3 iteration of the CMS HH Multilepton analysis.

The code is forked and for now heavily based on the UHH bersion of the [HH → bb𝜏𝜏 analysis](https://github.com/uhh-cms/hh2bbtautau)
and still very much WIP. Expect remnants from the bb𝜏𝜏 analysis, crashes and bugs, you have been warned!

Please make sure you are subscribed to our e-group: cms-hh-multilepton@cern.ch
It controls the acess to our indico etc. and is a good way to get updates for our meetings.

Also join our channel on [mattermost](https://mattermost.web.cern.ch/cms-exp/channels/hh-multilepton-run3).
(You will need to join the CMS team first if not done so).

The code is currently developed with the Tallinn T2 (and lxplus) in mind.
For further questions please, contact t\*\*\*\*.l\*\*\*\*@no-spam-cern.ch .

## First time setup

```shell
# clone the project
git clone --recursive git@github.com:HEP-KBFI/hhmultilepton.git
cd hhmultilepton

# source the setup and store decisions in .setups/dev.sh (arbitrary name)
source setup.sh dev

# Decisions include storage locations, these should be set according to the system you are running the code on:
# CF_DATA should point to a location in home (manivald) or afs (lxplus), same as CF_SOFTWARE_BASE and CF_JOB_BASE
# CF_WLCG_CACHE_ROOT is a cache for remote files should be on /local/user (manivald) or eos (lxplus).

# suggestion for lxplus setup:

export CF_DATA="$CF_REPO_BASE/data"
export CF_SOFTWARE_BASE="$CF_DATA/software"
export CF_JOB_BASE="$CF_DATA/jobs"
export CF_STORE_NAME="cf_store"
export CF_STORE_LOCAL="$CF_DATA/$CF_STORE_NAME"
export CF_WLCG_CACHE_ROOT="/eos/user/$CF_CERN_USER_FIRSTCHAR/$CF_CERN_USER/HHMultilepton_Run3/cf_scratch"
export CF_WLCG_USE_CACHE="true"
export CF_WLCG_CACHE_CLEANUP="false"
export CF_CRAB_STORAGE_ELEMENT="T2_EE_Estonia"
export CF_CRAB_BASE_DIRECTORY="/store/user/$CF_CERN_USER/cf_crab_outputs"

# if space in afs is a problem, one can try
export CF_SOFTWARE_BASE="/eos/user/$CF_CERN_USER_FIRSTCHAR/$CF_CERN_USER/HHMultilepton_Run3/software"

# suggestion for manivald is the same except:
export CF_JOB_BASE="/local/$CF_CERN_USER/HHMultilepton_Run3/jobs"
export CF_WLCG_CACHE_ROOT="/local/$CF_CERN_USER/HHMultilepton_Run3/cf_scratch"

# After first time setup, *if on manivald the estonian login node*, open the created setup file and add:
export TMPDIR="/scratch/local/$CF_CERN_USER"

# get a voms token:

voms-proxy-init -voms cms -rfc -valid 196:00
```

<img width="1336" height="506" alt="image" src="https://github.com/user-attachments/assets/29e6f810-e273-4b2e-9a80-02427e228298" />


Code can now be run but first storage locations for the tasks outputs should be checked as configured [here](https://github.com/HEP-KBFI/hhmultilepton/blob/master/law_outputs.cfg#L26-L90)
Currently outputs point to the user store of the T2 on manivald so that outputs are also accessible remotely, but we will likely adapt this over time depending on the output.
I.e large outputs available in a remote reachable location, smaller ones on local stores. Larger ones likely also split by user/cluster so that central versions can be reused.

*Important* For development on lxplus *i strongly * advise to change wlcg_fs_manivald to wlcg_fs_cernbox in the beginning.

After this is set, try to run on signal locally:

```shell
law run cf.PlotVariables1D \
    --version test \
    --producers default \
    --variables nmu \
    --datasets hh_ggf_htt_hvv_kl1_kt1_powheg \
```

And if this runs on background via slurm/condor

```shell
law run cf.PlotVariables1D \
    --version test \
    --producers default \
    --variables nmu \
    --datasets zz_pythia \
    --workflow slurm \
```

or with

```shell
    --workflow htcondor \
```

crab to be tested.

## Documentation
TODO but a general overview can be found in these slides: https://indico.cern.ch/event/1580193/contributions/6660044/attachments/3121091/5534653/multilep%20framework.pdf

## 🙏 Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/tolange"><img src="https://avatars.githubusercontent.com/u/11850680?s=96&v=4" width="100px;" alt="`Torben Lange`"/><br /><sub><b>Torben Lange</b></sub></a><br /><a href="https://github.com/HEP-KBFI/hhmultilepton/commits/master/?author=tolange" title="Code">💻</a> </td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/MatheuspCoelho"><img src="https://avatars.githubusercontent.com/u/85200761?v=4" width="100px;" alt="`Matheus Coelho`"/><br /><sub><b>Matheus Coelho</b></sub></a><br /><a href="https://github.com/HEP-KBFI/hhmultilepton/commits/master/?author=MatheuspCoelho" title="Code">💻</a> </td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->


## Useful links

- [columnflow documentation](https://columnflow.readthedocs.io/en/latest/index.html)
- CMS services
  - [HLT info browser](https://cmshltinfo.app.cern.ch/path/HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v)
  - [HLT config browser](https://cmshltcfg.app.cern.ch/open?db=online&cfg=%2Fcdaq%2Fphysics%2FRun2018%2F2e34%2Fv2.1.5%2FHLT%2FV2)
  - [GrASP](https://cms-pdmv-prod.web.cern.ch/grasp/)
  - [XSDB](https://xsdb-temp.app.cern.ch)
  - [DAS](https://cmsweb.cern.ch/das)
NanoAOD:
  - [Nano documentation](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc)
  - [Correctionlib files](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration)
- JME
  - [Docs](https://cms-jerc.web.cern.ch)
- BTV
  - [Docs](https://btv-wiki.docs.cern.ch)
- TAU
  - [Run 2 Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2)
  - [Run 3 Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun3)
  - [Correctionlib files](https://gitlab.cern.ch/cms-tau-pog/jsonpog-integration/-/tree/TauPOG_v2_deepTauV2p5/POG/TAU?ref_type=heads)

## Development

- Source hosted at [GitHub](https://github.com/HEP-KBFI/hhmultilepton)
- Report issues, questions, feature requests on [GitHub Issues](https://github.com/HEP-KBFI/hhmultilepton/issues)
- Ideally also ping us on MM
- For new features open a new branch before merging into master, ask for a code review by a felllow contributor and dont forget linting!
- Happy coding :)
