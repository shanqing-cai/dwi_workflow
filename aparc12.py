#!/usr/bin/env python

# 1st column: primary name of ROI
# 2nd column: cortical zone (lobe)
# 3rd column: N - non-speech; S - speech
# 4th column: alternative name(s), if any

aROIs = [   ['FP',     'Prefrontal',   'N', []], \
            ['aMFg',   'Prefrontal',   'N', []], \
            ['aIFs',   'Prefrontal',   'N', []], \
            ['SFg',    'Prefrontal',   'N', []], \
            ['aFO',    'Prefrontal',   'S', []], \
            ['FOC',    'Prefrontal',   'N', []], \
            ['pMFg',   'Prefrontal',   'N', []], \
            ['pIFs',   'Prefrontal',   'S', []], \
            ['dIFt',   'Prefrontal',   'N', []], \
            ['dIFo',   'Prefrontal',   'S', []], \
            ['vIFt',   'Prefrontal',   'N', []], \
            ['vIFo',   'Prefrontal',   'S', []], \
            ['pFO',    'Prefrontal',   'S', []], \
            ['FMC',    'Prefrontal',   'N', []], \
            ['aINS',   'Insular',      'S', []], \
            ['pINS',   'Insular',      'S', []], \
            ['adPMC',  'Premotor',     'N', []], \
            ['mdPMC',  'Premotor',     'N', []], \
            ['pdPMC',  'Premotor',     'S', []], \
            ['midPMC', 'Premotor',     'S', []], \
            ['vPMC',   'Premotor',     'S', []], \
            ['SMA',    'Premotor',     'S', []], \
            ['preSMA', 'Premotor',     'S', []], \
            ['aCO',    'Precentral',   'S', []], \
            ['dMC',    'Precentral',   'S', []], \
            ['midMC',  'Precentral',   'S', []], \
            ['vMC',    'Precentral',   'S', []], \
            ['pCO',    'Postcentral',  'S', []], \
            ['dSC',    'Postcentral',  'S', []], \
            ['vSC',    'Postcentral',  'S', []], \
            ['aSMg',   'PPC',          'S', []], \
            ['pSMg',   'PPC',          'N', []], \
            ['PO',     'PPC',          'S', []], \
            ['SPL',    'PPC',          'N', []], \
            ['Ag',     'PPC',          'N', ["AG"]], \
            ['PCN',    'PPC',          'N', []], \
            ['TP',     'Temporal',     'N', []], \
            ['PP',     'Temporal',     'S', []], \
            ['Hg',     'Temporal',     'S', ["H"]], \
            ['PT',     'Temporal',     'S', []], \
            ['aSTg',   'Temporal',     'S', []], \
            ['pSTg',   'Temporal',     'S', []], \
            ['pdSTs',  'Temporal',     'S', []], \
            ['adSTs',  'Temporal',     'N', []], \
            ['pvSTs',  'Temporal',     'N', []], \
            ['avSTs',  'Temporal',     'N', []], \
            ['aMTg',   'Temporal',     'N', []], \
            ['pMTg',   'Temporal',     'N', []], \
            ['pITg',   'Temporal',     'N', []], \
            ['aCGg',    'Cingulate',    'S', ["aCG"]], \
            ['midCGg',  'Cingulate',    'N', ["midCG"]], \
            ['pCGg',    'Cingulate',    'N', ["pCG"]], \
            ['OC',     'Occipital',     'N', []], \
            ['MTOg',    'Occipital',    'N', ["MTO"]], \
            ['ITOg',    'Occipital',    'N', ["ITO"]], \
            ['Lg',     'Occipital',     'N', ["LG"]], \
            ['pPHg',    'Parahippocampal',  'N', ["pPH"]], \
            ['aPHg',    'Parahippocampal',  'N', ["aPH"]], \
            ['SCC',    'Parahippocampal',   'N', []]]
