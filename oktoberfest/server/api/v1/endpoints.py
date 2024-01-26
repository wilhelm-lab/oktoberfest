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
    m = sha256()
    m.update(hashInput.encode("utf-8"))
    m.update(SALT.encode("utf-8"))
    return m.hexdigest()[:HASHSIZE]


def dbGetStatus(taskId):
    with sqlite3.connect(DB_FILE) as conn:
        c = conn.cursor()
        c.execute("SELECT STATUS FROM JOBS WHERE TASK_ID = ?", [taskId])
        return c.fetchone()


def dbFilesUploaded(taskId):
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
    return {"status": 200}


@v1_blueprint.route("/submitJob", methods=["POST"])
def submitJob():
    taskId = calculate_hash(request.args.get("hashInput"))
    jobConfig = request.json
    jobConfig["taskId"] = taskId
    try:
        dbSubmitJob(taskId, json.dumps(jobConfig))
    except TypeError:
        return {"message": "You got a hash collision. You should increase hashsize."}
    filePath = os.path.join(WORKDIR, taskId, "config.json")
    with open(filePath, "w") as f:
        f.write(json.dumps(jobConfig))
    return jobConfig


@v1_blueprint.route("/uploadFile", methods=["POST"])
def uploadFile():
    taskId = calculate_hash(request.args.get("hashInput"))
    folderPath = os.path.join(WORKDIR, taskId)

    if request.args.get("taskId") != taskId:
        # Only for first file uploaded
        try:
            os.mkdir(folderPath)
            dbFilesUploaded(taskId)
        except FileExistsError:
            return {"message": "You got a hash collision. You should increase hashsize."}

    file = request.files["file"]
    filePath = os.path.join(folderPath, secure_filename(file.filename))
    file.save(filePath)

    return {"taskId": taskId}


@v1_blueprint.route("/jobStatus")
def jobStatus():
    taskId = calculate_hash(request.args.get("hashInput"))
    resp = {"id": taskId}
    try:
        resp["status"] = dbGetStatus(taskId)[0]
    except TypeError:
        resp["status"] = "UNKOWN"
    return resp


@v1_blueprint.route("/downloadResults")
def downloadResults():
    taskId = calculate_hash(request.args.get("hashInput"))
    folder = os.path.abspath(os.path.join(WORKDIR, taskId))
    if os.path.isfile(os.path.join(folder, "results.zip")):
        return send_from_directory(folder, "results.zip", as_attachment=True)
    else:
        return {"message": "Results not yet ready"}


@v1_blueprint.route("/getModels")
def getModels():
    return {
        "intensity": [
            {"name": "Prosit_2019_intensity"},
            {"name": "Prosit_2020_intensity_CID"},
            {"name": "Prosit_2020_intensity_HCD"},
            {"name": "Prosit_2020_intensity_TMT"},
            {"name": "Prosit_2020_intensity_CID"},
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
