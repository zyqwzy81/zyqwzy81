#!/bin/bash
# 血磷论文文件一键上传Github脚本
# 仓库本地路径
LOCAL_PATH="/data/小米云盘/文档/论文/血磷/github_repo"
# Github仓库HTTPS地址
GIT_URL="https://github.com/zyqwzy81/zyqwzy81.git"
# 提交备注
COMMIT_MSG="上传血磷论文相关全部文件"
# 分支名（你的本地是master）
BRANCH="master"

echo "==================== Git上传脚本开始 ===================="
# 切换到目标目录
cd "$LOCAL_PATH" || { echo "错误：目录不存在！$LOCAL_PATH"; exit 1; }
echo "当前工作目录：$(pwd)"

# 1. 清理原有远程origin
echo "1. 清理旧远程配置..."
git remote remove origin 2>/dev/null

# 2. 初始化仓库
echo "2. 初始化Git仓库..."
git init

# 3. 绑定Github远程仓库
echo "3. 绑定Github远程地址：$GIT_URL"
git remote add origin "$GIT_URL"

# 4. 添加所有文件
echo "4. 添加目录内全部文件至暂存区..."
git add .

# 5. 提交文件
echo "5. 提交文件，备注：$COMMIT_MSG"
git commit -m "$COMMIT_MSG" || echo "提示：无新增文件，工作区干净"

# 6. 拉取远程代码避免冲突
echo "6. 拉取远程仓库代码..."
git pull origin "$BRANCH" --allow-unrelated-histories

# 7. 推送至Github
echo "7. 推送文件到Github $BRANCH 分支"
echo "⚠️ 弹窗输入账号：zyqwzy81，密码处粘贴Github Token"
git push origin "$BRANCH"

echo "==================== 脚本执行完成 ===================="
