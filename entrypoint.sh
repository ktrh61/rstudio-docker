#!/bin/bash
set -e

# 環境変数からUID/GIDを取得（デフォルト値あり）
TARGET_UID=${HOST_UID:-$(id -u rstudio)}
TARGET_GID=${HOST_GID:-$(id -g rstudio)}

# 現在のrstudioユーザーのUID/GIDを確認
CURRENT_UID=$(id -u rstudio)
CURRENT_GID=$(id -g rstudio)

# UID/GIDが異なる場合のみ変更実行
if [ "$TARGET_UID" != "$CURRENT_UID" ] || [ "$TARGET_GID" != "$CURRENT_GID" ]; then
    echo "Adjusting UID/GID: $CURRENT_UID:$CURRENT_GID -> $TARGET_UID:$TARGET_GID"
    
    # グループ変更
    if [ "$TARGET_GID" != "$CURRENT_GID" ]; then
        groupmod -g "$TARGET_GID" rstudio
    fi
    
    # ユーザー変更
    if [ "$TARGET_UID" != "$CURRENT_UID" ]; then
        usermod -u "$TARGET_UID" rstudio
    fi
    
    # ホームディレクトリの所有権修正
    chown -R rstudio:rstudio /home/rstudio
    
    echo "UID/GID adjustment completed"
else
    echo "UID/GID already correct: $TARGET_UID:$TARGET_GID"
fi

# 元のエントリーポイント実行
exec /init "$@"
