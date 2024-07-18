"""
conda activate esmfold
pip install flask
pip install gunicorn
gunicorn -w 1 -b 127.0.0.1:8000 local_api:app
"""
from flask import Flask, request, jsonify
from evodesign.Prediction.ESMFold import ESMFold
from evodesign.Metrics.ESM2Descriptors import ESM2Descriptors
import json



esmfold_model = ESMFold()
# cargamos el modelo a la GPU de una vez
esmfold_model.predict_pdb_str('GREETINGS')
esm_model = ESM2Descriptors()
# cargamos el modelo a la GPU de una vez
esm_model.compute_descriptor_vectors('GREETINGS')





app = Flask(__name__)

@app.route('/esmfold', methods=[ 'POST' ])
def esmfold():
    data = request.get_json()
    sequence = data['sequence']
    prediction = esmfold_model.predict_pdb_str(sequence)
    return jsonify({ 'pdb': prediction })



@app.route('/esm', methods=[ 'POST' ])
def esm():
    data = request.get_json()  # Get data from the request
    sequence = data['sequence']
    matrix = esm_model.compute_descriptor_vectors(sequence)
    return jsonify({ 'descriptors': json.dumps(matrix.tolist()) })



if __name__ == '__main__':
    app.run(debug=False)
