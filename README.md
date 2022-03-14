# COMAD: COMmunity Assembly Dynamics of microbes

## Installation 

```bash
pip install -e.
```

## Sample run starting with biom file
```bash
comad neufit_biom --biom comad/tests/data/sample_biom --output_filename github_example --output_filepath comad/tests/data/testing_output
```

## Sample run starting with data and taxonomy files
```bash
comad neufit --_data_filename comad/tests/data/sample_data.csv --_taxonomy_filename comad/tests/data/sample_taxonomy.csv --output_filename github_example --output_folder_path comad/tests/data/testing_output/github_example
```

Modified by: Caitlin Guccione
    
    Copyright
    ---------
    Github: https://github.com/misieber/neufit
    Theroy as described in [1]_. and [2]_.
    
    Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions 
    are met:

    1. Redistributions of source code must retain the above copyright 
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above 
       copyright notice, this list of conditions and the following 
       disclaimer in the documentation and/or other materials provided 
       with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
    
    
    References
    ----------
    .. [1] Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S. 
    and Curtis, T. P. (2006). Quantifying the roles of immigration 
    and chance in shaping prokaryote community structure. Environmental 
    Microbiology, 8:732-740. https://doi.org/10.1111/j.1462-2920.2005.00956.x
    
    .. [2] Sieber, M., Pita, L., Weiland-Bräuer, N., Dirksen, P., Wang, J., 
    Mortzfeld, B., Franzenburg, S., Schmitz, R. A., Baines, J. F., Fraune, 
    S., Hentschel, U., Schulenburg, H., Bosch, T. C. G. and Traulsen, A. 
    (2018). The Neutral Metaorganism. bioRxiv. https://doi.org/10.1101/367243
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
