"""
gunicorn -w 1 -b 127.0.0.1:8000 local_api_esm_only:app
"""

from flask import Flask, request, jsonify
from Metrics.ESM2 import ESM2


esm_model = ESM2(gpu_device="cuda:0")
esm_model.query_model("GREETINGS")


app = Flask(__name__)


@app.route("/esm", methods=["POST"])
def esm():
    data = request.get_json()
    sequence = data["sequence"]
    desc_matrix, predicted_contacts = esm_model.query_model(sequence)
    return jsonify(
        {
            "desc_matrix": desc_matrix.tolist(),
            "predicted_contacts": predicted_contacts.tolist(),
        }
    )


if __name__ == "__main__":
    app.run(debug=False)
