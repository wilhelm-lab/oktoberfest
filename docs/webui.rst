Web UI
======

Oktoberfest comes with a browser-based interface that lets you submit and monitor jobs through a guided form instead of writing JSON configuration files by hand.
The web UI runs locally on your machine (Windows, macOS and Linux) via Docker, or as a shared server for a group, and executes the same jobs as the high-level API (see `Running Oktoberfest <./jobs.html>`_).
Rescoring jobs include MS1 quantification (search engine or FlashLFQ) and PTM site localization with false localization rate estimation.
The Docker image contains Oktoberfest, Percolator, ThermoRawFileParser, FlashLFQ and the .NET runtime, so only the predictions are requested from a Koina server.

Supported jobs
--------------

.. table::
    :class: fixed-table

    +------------------------------+-----------------------------+---------------------------------------------------------------------+
    | Job type                     | Required inputs             | Outputs                                                             |
    +==============================+=============================+=====================================================================+
    | Rescoring                    | Search results + spectra    | Rescored PSM / peptide tables and summary plots                     |
    +------------------------------+-----------------------------+---------------------------------------------------------------------+
    | Collision Energy Calibration | Search results + spectra    | Calibrated collision energy and spectral angle plots                |
    +------------------------------+-----------------------------+---------------------------------------------------------------------+
    | Spectral Library Generation  | FASTA or peptides CSV       | Spectral library in msp, spectronaut or dlib format                 |
    +------------------------------+-----------------------------+---------------------------------------------------------------------+

.. note::
    Training data generation, local DLOmix models and FDR estimation with mokapot are currently not available in the web UI. Use the high-level API for these (see `Running Oktoberfest <./jobs.html>`_).

Running the web UI with Docker
------------------------------

Prerequisites:
  - `Docker <https://docs.docker.com/get-docker/>`_ (Docker Desktop on Windows / macOS, Docker Engine with Docker Compose on Linux). On Windows, WSL2 is recommended as the backend for Docker Desktop.
  - `Git <https://git-scm.com/downloads>`_

Clone the repository, change to the ``webui`` folder and create your settings file from the provided template:

.. code-block:: bash

   git clone https://github.com/wilhelm-lab/oktoberfest.git
   cd oktoberfest/webui
   cp oktoberfest-webui.env.example oktoberfest-webui.env
   mkdir -p data

Every setting is explained in ``oktoberfest-webui.env``; the defaults work for a personal installation. On Linux, set ``APP_UID`` and ``APP_GID`` to the output of ``id -u`` and ``id -g`` so that the job files belong to you.
By default, all uploaded inputs, generated outputs, the job database and logs are saved in ``webui/data/``. To store them somewhere else, e.g. on an external drive, change ``HOST_DATA_DIR``.

Build and start the containers in the background (always pass ``--env-file``, also for ``down`` and ``logs``):

.. code-block:: bash

   docker compose --env-file oktoberfest-webui.env up --build -d

Then open `http://localhost:8080 <http://localhost:8080>`_ in your browser. ``scripts/smoke_test.sh`` checks that everything works. Jobs and data are kept in ``HOST_DATA_DIR`` and survive container restarts. To stop the web UI, run:

.. code-block:: bash

   docker compose --env-file oktoberfest-webui.env down

Submitting and monitoring jobs
------------------------------

Select the job type, upload the inputs and fill in the form. After submitting, you are redirected to ``/jobs/<job_id>`` where the status is updated automatically. Large jobs can take hours, so bookmark this page and come back later. All recent jobs are listed at ``/jobs``.

Jobs are executed one after the other by default. Each job gets its own directory in ``DATA_DIR``:

.. code-block:: text

   <job_id>/
     inputs/        uploaded files
     output/        Oktoberfest output (see Generated outputs)
     config.json    the resolved Oktoberfest configuration
     captured.log   log of the Oktoberfest run
     results.zip    downloadable archive of the results and config

The ``config.json`` of a job can be reused to run the same job with the high-level API. The content of ``output/`` is described in `Generated outputs <./outputs.html>`_.

Configuration
-------------

The web UI is configured in ``oktoberfest-webui.env``. The most important settings:

.. table::
    :class: fixed-table

    +----------------------+------------------------------+------------------------------------------------------------------------+
    | Variable             | Default                      | Description                                                            |
    +======================+==============================+========================================================================+
    | HOST_DATA_DIR        | ./data                       | Folder that stores all jobs, uploads, results and the job database     |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | APP_UID / APP_GID    | 1000                         | User and group the containers run as (owner of the job files)          |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | API_PORT             | 8080                         | Port of the web UI                                                     |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | APP_MODE             | local                        | local (single user) or hosted (shared server with separate sessions)   |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | CELERY_CONCURRENCY   | 1                            | Number of jobs executed in parallel                                    |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | MAX_TOTAL_THREADS    | 4                            | CPU threads all running jobs may use together                          |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | PREDICTION_SERVER    | koina.wilhelmlab.org:443     | Default Koina server, can also be changed per job in the form          |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | MAX_UPLOAD_BYTES     | 5000000000                   | Maximum size per uploaded file (5 GB)                                  |
    +----------------------+------------------------------+------------------------------------------------------------------------+

.. note::
    Full proteome rescoring jobs can require 8-16 GB of RAM depending on the size of the dataset. Each job uses one thread per spectra file (at most 4).

    Running the web UI as a shared server (hosted mode) needs HTTPS and an IP salt file in the data folder; see the web UI README.

For local development without Docker, running the tests and details about the architecture, see the `web UI README <https://github.com/wilhelm-lab/oktoberfest/blob/development/webui/README.md>`_.
