# ============================================================
# GitHub 推送脚本 (PowerShell) - 将 github_repo 推送到 GitHub
# 使用方法：
#   1. 将 github_repo 文件夹复制到你的电脑
#   2. 打开 PowerShell 进入 github_repo 目录
#   3. 运行: .\push_to_github.ps1
# ============================================================

$GITHUB_USER = "zyqwzy81"
$REPO_NAME = "zyqwzy81"
$REMOTE_SSH = "git@github.com:${GITHUB_USER}/${REPO_NAME}.git"
$REMOTE_HTTPS = "https://github.com/${GITHUB_USER}/${REPO_NAME}.git"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host " GitHub 推送脚本" -ForegroundColor Cyan
Write-Host " 仓库: ${GITHUB_USER}/${REPO_NAME}" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan

# 检查是否在 git 仓库内
if (-not (Test-Path ".git")) {
    Write-Host "[错误] 当前目录不是 Git 仓库" -ForegroundColor Red
    Write-Host "请确保在 github_repo 目录下运行此脚本" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "[1/3] 设置远程仓库..." -ForegroundColor Green

# 默认使用 SSH（如果已经配置了 SSH 密钥）
git remote set-url origin $REMOTE_SSH

Write-Host ""
Write-Host "[2/3] 推送代码到 GitHub..." -ForegroundColor Green
Write-Host "     目标: $REMOTE_SSH" -ForegroundColor Gray

try {
    git push -u origin master 2>&1 | Out-Null
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  -> 推送成功！" -ForegroundColor Green
    } else {
        throw "Push failed"
    }
} catch {
    Write-Host ""
    Write-Host "[!] SSH 推送失败，尝试 HTTPS..." -ForegroundColor Yellow
    git remote set-url origin $REMOTE_HTTPS
    
    Write-Host "正在使用 HTTPS 推送（需要输入认证信息）..." -ForegroundColor Yellow
    Write-Host "密码/Token 输入方法：" -ForegroundColor Yellow
    Write-Host "  1. 在 https://github.com/settings/tokens 生成 Personal Access Token" -ForegroundColor Yellow
    Write-Host "  2. 勾选 repo 权限，复制 token" -ForegroundColor Yellow
    Write-Host "  3. 在弹出窗口的用户名输入 ${GITHUB_USER}，密码输入 token" -ForegroundColor Yellow
    Write-Host ""
    
    git push -u origin master
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  -> 推送成功！" -ForegroundColor Green
    } else {
        Write-Host "[错误] 推送失败" -ForegroundColor Red
        Write-Host ""
        Write-Host "手动推送命令：" -ForegroundColor Yellow
        Write-Host "  git remote set-url origin $REMOTE_SSH" -ForegroundColor Gray
        Write-Host "  git push -u origin master" -ForegroundColor Gray
        exit 1
    }
}

Write-Host ""
Write-Host "[3/3] 完成！" -ForegroundColor Green
Write-Host "  -> 打开 https://github.com/${GITHUB_USER}/${REPO_NAME} 查看" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan