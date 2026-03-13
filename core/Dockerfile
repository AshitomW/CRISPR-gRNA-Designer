FROM python:3.14-alpine AS builder
WORKDIR /app
RUN apk add --no-cache curl build-base
RUN pip install --no-cache-dir poetry
RUN poetry config virtualenvs.create false
COPY pyproject.toml poetry.lock* ./
RUN poetry install --no-interaction --no-ansi --no-root
COPY src/ src/

FROM python:3.14-alpine
WORKDIR /app
RUN apk add --no-cache curl
COPY --from=builder /usr/local/lib/python3.14/site-packages /usr/local/lib/python3.14/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin
COPY --from=builder /app/src/ src/
ENV PYTHONPATH=/app/src
CMD ["poetry", "run", "uvicorn", "crispr_designer.main:app", "--host", "0.0.0.0", "--port", "8000"]
