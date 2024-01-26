from flask import Blueprint, request, session, send_from_directory
from hashlib import sha256
import datetime
from ...config import WORKDIR, HASHSIZE, SALT, DB_FILE
from werkzeug.utils import secure_filename
import os
import json

import sqlite3


v1_blueprint = Blueprint("v1", __name__)


def calculate_hash(hashInput):
    """
    Calculate a SHA256 hash of the input string, salted with a global salt value.

    This function takes an input string, encodes it to bytes using UTF-8 encoding, and then 
    calculates a SHA256 hash of the encoded string. The global salt value is also encoded to 
    bytes and added to the hash. The resulting hash is then truncated to the size specified 
    by the global HASHSIZE value.

    :param hashInput: The input string to hash.
    :type hashInput: str
    :return: The calculated hash.
    :rtype: str
    """
    m = sha256()
    m.update(hashInput.encode("utf-8"))
    m.update(SALT.encode("utf-8"))
    return m.hexdigest()[:HASHSIZE]


def dbGetStatus(taskId):
    """
    Get the status of a task from the database.

    This function connects to the database and retrieves the status of the job with the given taskId 
    from the JOBS table. If no status is found, None is returned.

    :param taskId: The ID of the task.
    :type taskId: str
    :return: The status of the task, or None if no status is found.
    :rtype: str or None
    """
    with sqlite3.connect(DB_FILE) as conn:
        c = conn.cursor()
        c.execute("SELECT STATUS FROM JOBS WHERE TASK_ID = ?", [taskId])
        return c.fetchone()


def dbFilesUploaded(taskId):
    """
    Update the database to indicate that files have been uploaded for a task.

    This function connects to the database, starts an exclusive transaction, and inserts a new 
    record into the JOBS table with the status 'UPLOADED'. The ID of the new record is one more 
    than the maximum ID currently in the JOBS table. The current date and time is also recorded. 
    If the database update is successful, the transaction is committed.

    :param taskId: The ID of the task.
    :type taskId: str
    """
    with sqlite3.connect(DB_FILE) as conn:
        c = conn.cursor()
        c.execute("BEGIN EXCLUSIVE")
        c.execute("SELECT MAX(ID) FROM JOBS")
        maxId = c.fetchone()[0]
        c.execute(
            "INSERT INTO JOBS (ID, TASK_ID, DATE, STATUS) VALUES (?, ?, ?, 'UPLOADED')",
            [maxId + 1, taskId, datetime.datetime.now()],
        )
        c.execute("COMMIT")


def dbSubmitJob(taskId, config):
    """
    Update the database to indicate that a job has been submitted and is ready to be executed.

    This function connects to the database, starts an exclusive transaction, and updates the 
    status of the job with the given taskId to 'PENDING'. It also inserts the job configuration 
    into the CONFIGS table. If the job submission is successful, the transaction is committed. 
    If the job submission fails, the transaction is rolled back.

    :param taskId: The ID of the task.
    :type taskId: str
    :param config: The configuration for the job.
    :type config: str
    """
    with sqlite3.connect(DB_FILE) as conn:
        c = conn.cursor()
        c.execute("BEGIN EXCLUSIVE")
        c.execute("SELECT ID FROM JOBS WHERE TASK_ID = ? AND STATUS = 'UPLOADED'", [taskId])
        dbId = c.fetchone()[0]
        c.execute("UPDATE JOBS SET STATUS = 'PENDING' WHERE ID = ?", [dbId])
        c.execute("INSERT INTO CONFIGS (ID, CONFIG) VALUES (?, ?)", [dbId, config])
        c.execute("COMMIT")


@v1_blueprint.route("/status")
def status():
    """
    Endpoint to check if the server is running.
    """
    return "Running", 200

@v1_blueprint.route("/runningJobs")
def runningJobs():
    return {"count": 5}

@v1_blueprint.route("/submitJob", methods=["POST"])
def submitJob():
    """
    Endpoint to submit a job.

    This endpoint accepts a JSON object in the request body containing the job configuration. 
    The 'hashInput' parameter in the request args is used to calculate the taskId. The job 
    configuration and taskId are then used to submit the job to the database. If the job 
    submission is successful, the config written to the database is returned. If the job submission fails, 
    an error message is returned.

    :param hashInput: The input string used to calculate the taskId.
    :type hashInput: str
    :param config: The configuration for the job, passed in the request body.
    :type config: dict
    :return: A dictionary containing a success or error message.
    :rtype: dict
    """
    hashInput = request.args.get("hashInput")
    if hashInput:
        taskId = calculate_hash(hashInput)
        jobConfig = request.json
        jobConfig["taskId"] = taskId
        try:
            dbSubmitJob(taskId, json.dumps(jobConfig))
        except TypeError:
            return "Hash collision detected. Reload and submit again.", 409
        filePath = os.path.join(WORKDIR, taskId, "config.json")
        with open(filePath, "w") as f:
            f.write(json.dumps(jobConfig))
        return jobConfig
    else:
        return "No hashInput provided", 400


@v1_blueprint.route("/uploadFile", methods=["POST"])
def uploadFile():
    """
    Endpoint to upload a file for a task. 

    This endpoint accepts a file in the request body and two parameters in the request args: 
    'hashInput' and 'taskId'. The 'hashInput' is used to calculate the taskId. If the 'taskId' 
    parameter does not match the calculated taskId, it is assumed that this is the first file 
    being uploaded for the task and a new directory is created for the task. If a directory 
    already exists with the same taskId, a hash collision is assumed and an error message is returned.

    :param hashInput: The input string used to calculate the taskId.
    :type hashInput: str
    :param taskId: The ID of the task. This should match the taskId calculated from the hashInput.
    :type taskId: str
    :param file: The file to be uploaded.
    :type file: werkzeug.datastructures.FileStorage
    :return: A dictionary containing the taskId.
    :rtype: dict
    """
    taskId = calculate_hash(request.args.get("hashInput"))
    folderPath = os.path.join(WORKDIR, taskId)

    if request.args.get("taskId") != taskId:
        # Only for first file uploaded
        try:
            os.mkdir(folderPath)
            dbFilesUploaded(taskId)
        except FileExistsError:
            return "Hash collision detected. Reload and submit again.", 409

    file = request.files["file"]
    filePath = os.path.join(folderPath, secure_filename(file.filename))
    file.save(filePath)

    return {"taskId": taskId}


@v1_blueprint.route("/jobStatus")
def jobStatus():
    """
    Endpoint to get the status of a job.

    This endpoint accepts a 'hashInput' parameter in the request args, which is used to calculate 
    the taskId. The status of the job with the calculated taskId is then retrieved from the database. 
    If no status is found, the status is set to 'UNKNOWN'.

    :param hashInput: The input string used to calculate the taskId.
    :type hashInput: str
    :return: A dictionary containing the taskId and the status of the job.
    :rtype: dict
    """
    taskId = request.args.get("taskId")
    if taskId:
        resp = {"id": taskId}
        try:
            resp["status"] = dbGetStatus(taskId)[0]
        except TypeError:
            resp["status"] = "UNKOWN"
        return resp
    else:
        return "No taskId provided", 400


@v1_blueprint.route("/downloadResults")
def downloadResults():
    """
    Endpoint to download the results of a job.

    This endpoint accepts a 'hashInput' parameter in the request args, which is used to calculate 
    the taskId. The results of the job with the calculated taskId are then retrieved from the 
    file system. If no results are found, an error message is returned.

    :param hashInput: The input string used to calculate the taskId.
    :type hashInput: str
    :return: A file containing the results of the job.
    :rtype: werkzeug.datastructures.FileStorage
    """
    taskId = request.args.get("taskId")
    if taskId:
        folder = os.path.abspath(os.path.join(WORKDIR, taskId))
        if os.path.isfile(os.path.join(folder, "results.zip")):
            return send_from_directory(folder, "results.zip", as_attachment=True)
        else:
            return "Result files not found", 404
    else:
        return "No taskId provided", 400


@v1_blueprint.route("/getModels")
def getModels():
    return {
        "intensity": [
            {"name": "Prosit_2019_intensity"},
            {"name": "Prosit_2020_intensity_CID"},
            {"name": "Prosit_2020_intensity_HCD"},
            {"name": "Prosit_2020_intensity_TMT"}
        ],
        "irt": [
            {"name": "Prosit_2019_irt"}, 
            {"name": "Prosit_2020_irt_TMT"}, 
            {"name": "Deeplc_hela_hf"}, 
            {"name": "AlphaPept_rt_generic"}
        ],
    }


@v1_blueprint.route("/getSearchEngines")
def getSearchEngines():
    return [{"name": "MaxQuant"}, {"name": "MSFragger"}, {"name": "Sage"}]


@v1_blueprint.route("/getRawFileTypes")
def getRawFileTypes():
    return [{"name": "thermo"}, {"name": "mzml"}, {"name": "bruker"}]
