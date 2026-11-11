FROM python:3.7.12-slim

WORKDIR /app

# Copy all project files
COPY . .

# pip automatically reads pyproject.toml, installs all dependencies, and installs evodesign
RUN pip install --no-cache-dir .

ENTRYPOINT ["python3", "-m", "evodesign"]
