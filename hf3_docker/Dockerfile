# Assume we installed the official AF3 Docker image
FROM alphafold3:latest

WORKDIR /app

# Install git if not already present
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# Clone HighFold 3
RUN git clone https://github.com/hongliangduan/HighFold3.git
WORKDIR /app/HighFold3
RUN git checkout 1e6125ad17993f4e951f577390da719c537cf01d
WORKDIR /app

# Overwrite AlphaFold 3 code with HighFold 3 modifications
# Adjust paths according to HighFold3 repo structure
RUN cp -r HighFold3/common/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/common
RUN cp -r HighFold3/constants/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/constants/
RUN cp -r HighFold3/data/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/data/
RUN cp -r HighFold3/jax/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/jax/
RUN cp -r HighFold3/model/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/model/
RUN cp -r HighFold3/parsers/cpp/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/parsers/cpp/
RUN cp -r HighFold3/scripts/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/scripts/
RUN cp -r HighFold3/structure/* /alphafold3_venv/lib/python3.12/site-packages/alphafold3/structure/
RUN cp -r HighFold3/run_alphafold.py /app/alphafold/run_alphafold.py

RUN pip install --no-cache-dir chex==0.1.88

WORKDIR /app/alphafold
