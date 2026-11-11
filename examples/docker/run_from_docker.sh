docker run --rm \
    -v "$(pwd)"/../:/data \
    -w /data \
    evodesign:latest \
    1y32.pdb \
    docker/settings.json \
    . \
    --jobname=test_docker
