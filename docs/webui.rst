Web UI
======

Oktoberfest comes with a browser-based interface that lets you submit and monitor jobs through a guided form instead of writing JSON configuration files by hand.
The web UI runs locally on your machine (Windows, macOS and Linux) via Docker and executes the same jobs as the high-level API (see `Running Oktoberfest <./jobs.html>`_).

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

Clone the repository, change to the ``webui`` folder and create your local environment file from the provided example:

.. code-block:: bash

   git clone https://github.com/wilhelm-lab/oktoberfest.git
   cd oktoberfest/webui
   cp .env.example .env

By default, all uploaded inputs, generated outputs, the job database and logs are saved in ``webui/data/``. To store them somewhere else, e.g. on an external drive, change ``DATA_DIR`` in the ``.env`` file.

Build and start the containers in the background:

.. code-block:: bash

   docker compose up --build -d

Then open `http://localhost:8000 <http://localhost:8000>`_ in your browser. Jobs and data are kept in ``DATA_DIR`` and survive container restarts. To stop the web UI, run:

.. code-block:: bash

   docker compose down

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

The web UI is configured via environment variables in the ``.env`` file:

.. table::
    :class: fixed-table

    +----------------------+------------------------------+------------------------------------------------------------------------+
    | Variable             | Default                      | Description                                                            |
    +======================+==============================+========================================================================+
    | DATA_DIR             | ./data                       | Root directory for job inputs, outputs and logs                        |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | DATABASE_URL         | sqlite:///./data/app.db      | Job database, SQLite (local) or Postgres URL                           |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | REDIS_URL            | redis://localhost:6379/0     | Message broker for the job queue                                       |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | PREDICTION_SERVER    | koina.wilhelmlab.org:443     | Default Koina server, can also be changed per job in the form          |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | MAX_UPLOAD_BYTES     | 5000000000                   | Maximum size per uploaded file (5 GB)                                  |
    +----------------------+------------------------------+------------------------------------------------------------------------+
    | CELERY_CONCURRENCY   | 1                            | Number of jobs executed in parallel (1 = one job at a time)            |
    +----------------------+------------------------------+------------------------------------------------------------------------+

.. note::
    ThermoRawFileParser is not included in the Docker image. Either provide the path to ThermoRawFileParser in the advanced section of the form, or convert your RAW files to mzML beforehand (e.g. with msConvert) and select mzML as spectra type.

    Full proteome rescoring jobs can require 8-16 GB of RAM depending on the size of the dataset. The number of processes used per job can be set via ``numThreads`` in the advanced section of the form.

For local development without Docker, running the tests and details about the architecture, see the `web UI README <https://github.com/wilhelm-lab/oktoberfest/blob/development/webui/README.md>`_.
