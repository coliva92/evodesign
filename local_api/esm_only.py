from flask import Flask, request, jsonify
from evodesign.Metrics.ESM2Descriptors import ESM2Descriptors



esm_model = ESM2Descriptors(gpu_device="cuda:0")
esm_model.compute_descriptor_vectors('GREETINGS')





app = Flask(__name__)

@app.route("/esm", methods=[ "POST" ])
def esm():
    data = request.get_json()
    sequence = data["sequence"]
    matrix = esm_model.compute_descriptor_vectors(sequence)
    return jsonify({ "descriptors": matrix.tolist() })



if __name__ == '__main__':
    app.run(debug=False)
