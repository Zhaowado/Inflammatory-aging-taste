####Figure3B####
y_prob = rf_final.predict_proba(X_test_scaled)

fpr = dict()
tpr = dict()
roc_auc = dict()
n_classes = len(classes)

y_test_bin = label_binarize(y_test, classes=classes)

for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_prob[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

plt.figure(figsize=(10, 8))
for i, class_name in enumerate(classes):
    plt.plot(fpr[i], tpr[i], lw=2, label=f'{class_name} (AUC = {roc_auc[i]:.3f})')

plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for Each Class')
plt.legend(loc="lower right")
plt.savefig('roc_curve.pdf', dpi=300)
plt.close()

####Figure3C####
cm = confusion_matrix(y_test, y_pred)
cm_percent = (cm / cm.sum(axis=1, keepdims=True) * 100)  # 每行百分比

plt.figure(figsize=(10, 8))
classes = ['P0', 'P14', 'M2', 'M12', 'M24']

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(
    cm_percent,
    annot=True,
    fmt=".1f",
    cmap=plt.cm.Blues,
    xticklabels=classes,
    yticklabels=classes,
    annot_kws={"fontsize": 14},
    ax=ax,
    cbar=True,
)
