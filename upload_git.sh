#!/bin/bash
# 血磷论文一键上传GitHub脚本
# -------------------------- 配置区（仅修改这里） --------------------------
GIT_USER="zyqwzy81"
LOCAL_PATH="/data/小米云盘/文档/论文/血磷/github_repo"
COMMIT_MSG="上传血磷论文相关全部文件"
BRANCH="master"
# -----------------------------------------------------------------------------

echo "==================== Git上传脚本开始 ===================="
# 进入工作目录，目录不存在直接退出
cd "$LOCAL_PATH" || { echo "错误：目录不存在 $LOCAL_PATH"; exit 1; }
echo "当前目录：$(pwd)"

# 优先使用 SSH 方式（无需 token，已在 ~/.ssh/config 配置 ssh.github.com:443）
# 若 origin 已存在，复用现有配置
if git remote get-url origin >/dev/null 2>&1; then
  echo "已存在 origin，复用现有远程配置"
else
  if [ -f ~/.ssh/id_ed25519.pub ]; then
    git remote add origin "git@github.com:${GIT_USER}/zyqwzy81.git"
    echo "已绑定远程仓库（SSH）: github.com/zyqwzy81/zyqwzy81"
  else
    # 备选：HTTPS + 环境变量 token（禁止把 token 写入脚本文件）
    if [ -z "$GITHUB_TOKEN" ]; then
      echo "错误：无 SSH key 且未设置环境变量 GITHUB_TOKEN（不要将 token 写入脚本）"
      echo "用法：export GITHUB_TOKEN=ghp_xxx ; ./upload_git.sh"
      exit 1
    fi
    git remote add origin "https://${GIT_USER}:${GITHUB_TOKEN}@github.com/${GIT_USER}/zyqwzy81.git"
    echo "已绑定远程仓库（HTTPS）: github.com/zyqwzy81/zyqwzy81"
  fi
fi

# 添加全部文件
git add .

# 提交代码，无文件提交时静默提示
git commit -m "$COMMIT_MSG" || echo "提示：工作区无新增修改文件"

# 拉取远程代码，合并无关历史
git pull origin "$BRANCH" --allow-unrelated-histories 2>/dev/null || true

# 推送至远程（如遇SSL错误，可取消下行注释）
# export GIT_SSL_NO_VERIFY=1
git push origin "$BRANCH"

echo "==================== 脚本执行结束 ===================="