#!/bin/bash
# ============================================================
# GitHub 推送脚本 - 将 github_repo 推送到 GitHub
# 使用方法：
#   1. 将 github_repo 文件夹复制到你的电脑
#   2. 打开终端进入 github_repo 目录
#   3. 运行: bash push_to_github.sh
# ============================================================

set -e

GITHUB_USER="zyqwzy81"
REPO_NAME="zyqwzy81"
REMOTE_SSH="git@github.com:${GITHUB_USER}/${REPO_NAME}.git"
REMOTE_HTTPS="https://github.com/${GITHUB_USER}/${REPO_NAME}.git"

echo "=========================================="
echo " GitHub 推送脚本"
echo " 仓库: ${GITHUB_USER}/${REPO_NAME}"
echo "=========================================="

# 检查是否在 git 仓库内
if [ ! -d ".git" ]; then
    echo "[错误] 当前目录不是 Git 仓库"
    echo "请确保在 github_repo 目录下运行此脚本"
    exit 1
fi

# 检查远程仓库是否可访问
echo ""
echo "[1/3] 检查远程仓库连接..."

# 优先尝试 SSH
if ssh -T -o StrictHostKeyChecking=no -o ConnectTimeout=5 "git@github.com" 2>&1 | grep -q "successfully"; then
    echo "  -> SSH 连接成功，使用 SSH 协议"
    git remote set-url origin "$REMOTE_SSH"
# 否则尝试 HTTPS
elif curl -s -o /dev/null -w "%{http_code}" --connect-timeout 5 "https://github.com" 2>/dev/null | grep -q "200"; then
    echo "  -> HTTPS 连接成功，使用 HTTPS 协议"
    git remote set-url origin "$REMOTE_HTTPS"
else
    echo "  -> 无法连接 GitHub，请检查网络"
    echo "  -> 将使用 HTTPS 协议（推送时可能需要输入用户名和 Personal Access Token）"
    git remote set-url origin "$REMOTE_HTTPS"
fi

echo ""
echo "[2/3] 推送代码到 GitHub..."
if git push -u origin master; then
    echo "  -> 推送成功！"
else
    echo ""
    echo "[!] 推送失败，常见原因："
    echo "  1. 没有 SSH 密钥 → 请先生成密钥并添加到 GitHub"
    echo "     ssh-keygen -t ed25519 -C \"your_email@example.com\""
    echo "     然后将 ~/.ssh/id_ed25519.pub 内容添加到："
    echo "     https://github.com/settings/ssh/new"
    echo ""
    echo "  2. 需要 Personal Access Token（HTTPS）→ 使用以下格式："
    echo "     用户名: ${GITHUB_USER}"
    echo "     密码:   在 https://github.com/settings/tokens 生成（勾选 repo 权限）"
    echo ""
    echo "  推送命令（手动重试）："
    echo "    git push -u origin master"
    exit 1
fi

echo ""
echo "[3/3] 验证推送结果..."
echo "  -> 打开 https://github.com/${GITHUB_USER}/${REPO_NAME} 查看"
echo ""
echo "=========================================="
echo " 完成！"
echo "=========================================="