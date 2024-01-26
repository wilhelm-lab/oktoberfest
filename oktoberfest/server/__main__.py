from flask import Flask, send_from_directory, abort, make_response, jsonify
from flask_cors import CORS, cross_origin
from .api.v1.endpoints import v1_blueprint
import time
from .config import MAX_CONTENT_LENGTH, UI_BUILD_DIR

app = Flask(__name__, static_folder='./dist', static_url_path='/prosit')
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH

app.register_blueprint(v1_blueprint, url_prefix='/api/v1')

@app.route('/')
def serve_vue_app():
    print(app.static_folder)
    return send_from_directory(app.static_folder, 'index.html')

if __name__ == '__main__':
    CORS(app, supports_credentials=True)
    app.config['CORS_HEADERS'] = 'Content-Type'
    app.run(debug=True)