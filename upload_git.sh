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

# 清除旧远程配置
git remote remove origin 2>/dev/null
# 初始化仓库
git init

# 使用国内镜像代理，解决TLS连接中断
AUTH_REMOTE="https://${GIT_USER}:${GIT_TOKEN}@github.com.cnpmjs.org/zyqwzy81/zyqwzy81.git"
git remote add origin "$AUTH_REMOTE"
echo "已绑定镜像代理远程仓库"


# 添加全部文件
git add .

# 提交代码，无文件提交时静默提示
git commit -m "$COMMIT_MSG" || echo "提示：工作区无新增修改文件"

# 拉取远程代码，合并无关历史
git pull origin "$BRANCH" --allow-unrelated-histories

# 推送至远程，无需手动输入账号密码
git push origin "$BRANCH"

echo "==================== 脚本执行结束 ===================="
