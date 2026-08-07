# RStudio Bioconductor Docker Makefile
# 複数環境対応（macOS + Colima、WSL2）

# 設定変数（既存run.shとの互換性保持）
IMG := bioconductor-openblas:3.21
PORT := 8787
PASSWORD := password
NAME := rstudio-cdm
THREADS := 1

# 環境変数での上書き対応
IMG := $(or $(IMG),bioconductor-openblas:devel)
PORT := $(or $(PORT),8787)
PASSWORD := $(or $(PASSWORD),password)
NAME := $(or $(NAME),rstudio-cdm)
THREADS := $(or $(THREADS),1)

# ユーザー情報（自動検出）
USER_UID := $(shell id -u)
USER_GID := $(shell id -g)
CURRENT_USER := $(shell whoami)

# 作業ディレクトリ
WORK_DIR := $(shell pwd)

# 環境検出とUID設定
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	OS_TYPE := macOS
	UID_ARGS := 
else
	OS_TYPE := Linux
	UID_ARGS := -e HOST_UID=$(USER_UID) -e HOST_GID=$(USER_GID)
endif

.PHONY: help build build-with-uid run stop clean logs info ps status

# デフォルトターゲット
help:
	@echo "=== RStudio Bioconductor Docker Makefile ==="
	@echo "Environment: $(OS_TYPE) | User: $(CURRENT_USER) ($(USER_UID):$(USER_GID))"
	@echo ""
	@echo "Available targets:"
	@echo "  build          - Build Docker image"
	@echo "  build-with-uid - Build with current user's UID/GID"
	@echo "  run            - Start RStudio container"
	@echo "  stop           - Stop RStudio container"
	@echo "  clean          - Stop and remove container"
	@echo "  logs           - Show container logs"
	@echo "  info           - Show session and system info"
	@echo "  ps             - Show container status"
	@echo ""
	@echo "Current settings:"
	@echo "  Image: $(IMG)"
	@echo "  Container: $(NAME)"
	@echo "  Port: $(PORT)"
	@echo "  Threads: $(THREADS)"

# 標準ビルド
build:
	@echo "Building Docker image: $(IMG)"
	docker build -t $(IMG) .
	@echo "✅ Build completed"

# UID指定ビルド（現在のユーザーに合わせる）
build-with-uid:
	@echo "Building Docker image with UID/GID: $(USER_UID)/$(USER_GID)"
	docker build \
		--build-arg USER_UID=$(USER_UID) \
		--build-arg USER_GID=$(USER_GID) \
		-t $(IMG) .
	@echo "✅ Build completed with custom UID/GID"

# コンテナ実行（既存run.shの機能を統合・macOS対応）
run:
	@echo "Starting RStudio container..."
	@echo "Environment: $(OS_TYPE) | User: $(CURRENT_USER)"
	@$(MAKE) stop 2>/dev/null || true
	docker run --rm -d \
		-e PASSWORD="$(PASSWORD)" \
		$(UID_ARGS) \
		-e OMP_NUM_THREADS="$(THREADS)" \
		-e OPENBLAS_NUM_THREADS="$(THREADS)" \
		-p "$(PORT):8787" \
		-v "$(WORK_DIR):/home/rstudio/project" \
		-w /home/rstudio/project \
		--name "$(NAME)" \
		"$(IMG)"
	@echo "✅ 起動: http://localhost:$(PORT)  user: rstudio  pass: $(PASSWORD)"

# コンテナ停止
stop:
	@echo "Stopping RStudio container..."
	@docker stop $(NAME) 2>/dev/null || echo "Container not running"

# クリーンアップ
clean: stop
	@echo "Removing RStudio container..."
	@docker rm $(NAME) 2>/dev/null || echo "Container already removed"

# ログ表示
logs:
	@echo "=== Container Logs ==="
	@docker logs $(NAME) 2>/dev/null || echo "Container not found"
	@echo ""
	@echo "=== Session Info ==="
	@grep -A3 -n "R version" out/logs/session_info.txt 2>/dev/null || echo "Session info file not found"

# システム情報
info:
	@echo "=== System Information ==="
	@echo "OS: $(OS_TYPE)"
	@echo "User: $(CURRENT_USER) (UID:$(USER_UID), GID:$(USER_GID))"
	@echo "Work Directory: $(WORK_DIR)"
	@echo "Docker Image: $(IMG)"
	@echo "Container Name: $(NAME)"
	@echo "Port: $(PORT)"
	@echo "Threads: $(THREADS)"
	@echo ""
	@echo "=== Docker Status ==="
	@docker --version 2>/dev/null || echo "Docker not available"
	@echo ""
	@echo "=== Container Status ==="
	@$(MAKE) ps

# コンテナ状態確認
ps:
	@docker ps -a --filter name=$(NAME) --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}" 2>/dev/null || echo "Docker not available"
