import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report, 
    confusion_matrix, 
    roc_curve, 
    auc, 
    roc_auc_score, 
    make_scorer
)
from sklearn.preprocessing import StandardScaler, label_binarize
# 设置随机种子确保结果可复现
np.random.seed(42)

adata = sc.read_h5ad('./00.data/taste.final.filter.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 去批次
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='timepoint')
sc.pp.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata)

# 识别DEG
sc.tl.rank_genes_groups(adata, groupby='timepoint', method='wilcoxon')

# 提取top基因
time_genes = []
for time in ['P0', 'P14', 'M2', 'M12', 'M24']:
    time_genes.extend(sc.get.rank_genes_groups_df(adata, group=time)['names'][:1000])

# 去除重复基因
time_genes = list(set(time_genes))

# 准备输入
X = adata[:, time_genes].X  # 基因表达矩阵
y = adata.obs['timepoint']  # 时间点标签

# 初始化五折分层交叉验证
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)

# 标准化
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 定义要评估的指标
scoring = {
    'accuracy': 'accuracy',
    'f1_macro': 'f1_macro',
    'roc_auc_ovr': make_scorer(
        roc_auc_score, 
        multi_class='ovr', 
        average='macro',
        needs_proba=True
    )
}

# 初始化随机森林模型
rf = RandomForestClassifier(n_estimators=1000, random_state=42, n_jobs=-1)

# 使用交叉验证评估模型
cv_results = cross_validate(
    rf, X_train_scaled, y_train, cv=cv, scoring=scoring, 
    return_train_score=True, 
    return_estimator=True,
    n_jobs=-1
)

# 输出交叉验证结果
print("五折交叉验证结果：")
print(f"训练集准确率: {np.mean(cv_results['train_accuracy']):.4f} ± {np.std(cv_results['train_accuracy']):.4f}")
print(f"测试集准确率: {np.mean(cv_results['test_accuracy']):.4f} ± {np.std(cv_results['test_accuracy']):.4f}")
print(f"训练集F1分数: {np.mean(cv_results['train_f1_macro']):.4f} ± {np.std(cv_results['train_f1_macro']):.4f}")
print(f"测试集F1分数: {np.mean(cv_results['test_f1_macro']):.4f} ± {np.std(cv_results['test_f1_macro']):.4f}")
print(f"测试集ROC AUC (OvR): {np.mean(cv_results['test_roc_auc_ovr']):.4f} ± {np.std(cv_results['test_roc_auc_ovr']):.4f}")

# 判断拟合状态
train_acc = np.mean(cv_results['train_accuracy'])
test_acc = np.mean(cv_results['test_accuracy'])
acc_diff = train_acc - test_acc

if train_acc < 0.7 and test_acc < 0.7:
    fit_status = "欠拟合"
elif acc_diff > 0.1:
    fit_status = "过拟合"
else:
    fit_status = "拟合良好"

print(f"模型状态: {fit_status} (训练-测试准确率差异: {acc_diff:.4f})")

# 使用全部训练数据重新训练最终模型
rf_final = RandomForestClassifier(n_estimators=1000, random_state=42, n_jobs=-1)
rf_final.fit(X_train_scaled, y_train)

# 模型评估（使用完整测试集）
y_pred = rf_final.predict(X_test_scaled)
print("\n最终模型在测试集上的表现：")
print(classification_report(y_test, y_pred))
