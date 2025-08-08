FROM bioconductor/bioconductor_docker:devel

# 環境変数（早い段階で設定）
ENV DEBIAN_FRONTEND=noninteractive \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# OpenBLAS とクリーンアップを一つのRUNで実行
RUN apt-get update && apt-get install -y --no-install-recommends \
      libopenblas0-pthread \
      libopenblas-dev \
      liblapack-dev \
      # 追加の有用なパッケージ
      htop \
      vim \
      curl \
      git \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get clean

# alternatives設定（エラーハンドリング強化）
RUN set -eux; \
    # アーキテクチャ検出
    ARCH=$(dpkg --print-architecture); \
    if [ "$ARCH" = "amd64" ]; then \
        ARCH_SUFFIX="x86_64-linux-gnu"; \
    elif [ "$ARCH" = "arm64" ]; then \
        ARCH_SUFFIX="aarch64-linux-gnu"; \
    else \
        echo "Unsupported architecture: $ARCH" >&2; \
        exit 1; \
    fi; \
    \
    # alternatives設定（アーキテクチャ対応）
    if update-alternatives --list libblas.so.3-${ARCH_SUFFIX} >/dev/null 2>&1; then \
        update-alternatives --set libblas.so.3-${ARCH_SUFFIX} /usr/lib/${ARCH_SUFFIX}/openblas-pthread/libblas.so.3; \
        update-alternatives --set liblapack.so.3-${ARCH_SUFFIX} /usr/lib/${ARCH_SUFFIX}/openblas-pthread/liblapack.so.3; \
    elif update-alternatives --list libblas.so.3 >/dev/null 2>&1; then \
        update-alternatives --set libblas.so.3 /usr/lib/${ARCH_SUFFIX}/openblas/libblas.so.3; \
        update-alternatives --set liblapack.so.3 /usr/lib/${ARCH_SUFFIX}/openblas/liblapack.so.3; \
    else \
        echo "Warning: Could not configure BLAS alternatives" >&2; \
    fi

# 動的UID/GID設定のための準備
ARG USER_UID=1000
ARG USER_GID=1000

# RStudioユーザーの権限設定（動的UID/GID対応）
RUN groupmod -g ${USER_GID} rstudio && \
    usermod -u ${USER_UID} rstudio && \
    chown -R rstudio:rstudio /home/rstudio

# エントリーポイントスクリプト追加
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# エントリーポイント設定
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

# 作業ディレクトリ設定
WORKDIR /home/rstudio

# ヘルスチェック追加
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s \
    CMD curl -f http://localhost:8787/ || exit 1

# RStudioサーバー起動（デフォルトコマンド）
CMD ["/init"]
