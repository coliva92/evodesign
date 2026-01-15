"""
gunicorn -w 1 -b 127.0.0.1:8000 local_api_esmfold_only:app
"""

from flask import Flask, request, jsonify
from evodesign.Prediction.ESMFold import ESMFold
from evodesign.Metrics.ESM2 import ESM2


gpu_device = "cuda:0"
esmfold_model = ESMFold(gpu_device)
esm_model = ESM2(gpu_device)
app = Flask(__name__)


@app.route("/esmfold", methods=["POST"])
def esmfold():
    data = request.get_json()
    sequence = data["sequence"]
    prediction = esmfold_model.predict_single_pdb_str(sequence)
    return jsonify({"pdb": prediction})


@app.route("/esm", methods=["POST"])
def esm_desc():
    data = request.get_json()
    sequence = data["sequence"]
    desc_matrix, predicted_contacts = esm_model.query_model(sequence)
    return jsonify(
        {
            "results": {
                "desc_matrix": desc_matrix.tolist(),
                "predicted_contacts": predicted_contacts.tolist(),
            }
        }
    )


if __name__ == "__main__":
    app.run(debug=False)
