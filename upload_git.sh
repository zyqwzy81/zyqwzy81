#!/bin/bash
# 血磷论文一键上传GitHub脚本
# -------------------------- 配置区（仅修改这里） --------------------------
GIT_USER="zyqwzy81"
# Personal access tokens (classic) 占位符，替换为你的真实ghp_开头token
GIT_TOKEN="YOUR_GITHUB_TOKEN"
LOCAL_PATH="/data/小米云盘/文档/论文/血磷/github_repo"
COMMIT_MSG="上传血磷论文相关全部文件"
BRANCH="master"
# -----------------------------------------------------------------------------

echo "==================== Git上传脚本开始 ===================="
# 进入工作目录，目录不存在直接退出
cd "$LOCAL_PATH" || { echo "错误：目录不存在 $LOCAL_PATH"; exit 1; }
echo "当前目录：$(pwd)"

# 若 origin 已存在（含已配置的凭据），复用现有配置，避免覆盖真实 token
if git remote get-url origin >/dev/null 2>&1; then
  echo "已存在 origin，复用现有远程配置（保留已配置凭据）"
else
  # 使用官方 GitHub 地址（如遇 TLS 连接问题，可取消注释下行 SSL 验证跳过）
  # 注意：github.com.cnpmjs.org 镜像代理已失效，不再使用
  AUTH_REMOTE="https://${GIT_USER}:${GIT_TOKEN}@github.com/zyqwzy81/zyqwzy81.git"
  git remote add origin "$AUTH_REMOTE"
  echo "已绑定远程仓库: github.com/zyqwzy81/zyqwzy81"
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