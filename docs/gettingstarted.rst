Getting started
**********************

Installation
============

Requirements:

* PySM `PySM <https://github.com/bthorne93/PySM_public>`_
* healpy

Clone the repository::

    pip install https://github.com/simonsobs/so_pysm_models/archive/master.zip

Development installation
========================

Clone from Github and install::

    git clone https://github.com/simonsobs/so_pysm_models
    cd so_pysm_models
    pip install -e .

Run unit tests::

    python setup.py test -V

Build docs::

    python setup.py build_docs -w
    
Example Usage
=============

This repository implements new models for PySM that can be added as additional components.

For example, create and configure a component::

    from so_pysm_models import GaussianSynchrotron
    synchrotron = GaussianSynchrotron(target_nside = 16)
    
Create a PySM sky and add this component::

    sky = pysm.Sky({})
    sky.add_component("gaussian_synch", gaussian_synch)

Then get a map at a specific frequency in GHz with standard PySM functionalities::

    m_synch = sky.gaussian_synch(2.3)

see example notebooks:

* `Example Gaussian Synchrotron <https://gist.github.com/zonca/51a6fa9763106c78813f964a4b88f0fc>`_
* `Example Gaussian Dust <https://gist.github.com/zonca/4ddb5e384cb34f8a2945c041d13e9428>`_
* `Example InterpolatingComponent <https://gist.github.com/zonca/08751497b040ec9d62ff5175573c786e>`_
* `Example Websky Extragalactic <https://gist.github.com/marcelo-alvarez/b13afc4a761c61334d11b5eeae953923>`_
