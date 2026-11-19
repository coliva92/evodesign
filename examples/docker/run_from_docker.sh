docker run --rm \
    -v "$(pwd)"/../:/data \
    -w /data \
    evodesign:latest \
    1y32.pdb \
    docker/settings.json \
    ../outputs \
    --jobname=test_docker
