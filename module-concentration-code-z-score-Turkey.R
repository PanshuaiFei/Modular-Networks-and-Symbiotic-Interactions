setwd("E:/Academic/SparCC/Gut_microbe_Protocol/新挑选微生物组22个/新挑选微生物组22个/data-for-run/total/network-cons/module含量")
rm(list = ls())
#读取 OTU 丰度表、样本分组和网络模块数据
otu<- read.csv('otu_divided.csv',check.names = F)

#otu=read.delim('otu.tsv',sep='\t',check.names=FALSE,row.names = 1)
#otu_divided <- otu / 10000000
# Save the data frame as a CSV file
#write.csv(otu_divided, file = "otu_divided.csv", row.names = TRUE)


group <- read.csv('group.csv')
module <- read.csv('node-module - 副本.csv')

#合并数据集
otu <- reshape2::melt(otu, id = 'otu')
names(otu) <- c('otu', 'sample', 'relative_abundance')
otu <- merge(otu, module, by = 'otu')
otu <- merge(otu, group, by = 'sample')
otu$module <- as.character(otu$module)
otu$group <- factor(otu$group)
head(otu)  #合并后的 OTU 在各样本的丰度、OTU 所属模块以及样本分组





library(dplyr)

# 删除 relative_abundance 为 0 的行
otu <- otu %>% filter(relative_abundance != 0)

# 如果您还想删除 NA 值，可以添加以下代码
otu <- otu %>% filter(!is.na(relative_abundance))

# 检查结果
head(otu)






#绘制绘制箱线图分别展示各模块各组微生物丰度的整体分布
library(ggplot2)

p <- ggplot(otu, aes(group, relative_abundance)) +
  geom_boxplot(aes(fill = module), show.legend = FALSE) +
  facet_wrap(~module, ncol = 3, scale = 'free_y') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'white')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = 'Group', y = 'Relative Abundance')


p

# 加载所需的包
library(ggplot2)
library(dplyr)
library(reshape2)

# 计算每个OTU在每个样本中的z-score
otu$z_score <- ave(otu$relative_abundance, otu$otu, FUN = function(x) (x - mean(x)) / sd(x))

# 对于每个模块，计算属于该模块的所有OTU的z-score的平均值
otu$module_z_score <- ave(otu$z_score, otu$sample, otu$module, FUN = mean)

# 使用ggplot2绘制箱线图
p2 <- ggplot(otu, aes(x = group, y = module_z_score)) +
  geom_boxplot(aes(fill = module), show.legend = FALSE) +
  facet_wrap(~module, ncol = 3, scale = 'free_y') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'white')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = 'Group', y = 'Relative Abundance Z-Score') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 显示图形
print(p2)











write.csv(otu,"z-score.csv")


p2 <- ggplot(otu, aes(x = group, y = module_z_score)) +
  geom_boxplot(aes(fill = module), show.legend = FALSE, outlier.shape = NA) +  # 添加outlier.shape = NA以隐藏离群点
  facet_wrap(~module, ncol = 3, scale = 'free_x') +  # 将scale设置为'free_x'，因为我们将手动设置y轴范围
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'white')) +
  scale_y_continuous(limits = c(-0.6, 2), expand = expansion(mult = c(0, 0))) +  # 明确设置y轴范围为0到2
  labs(x = 'Group', y = 'Relative Abundance Z-Score') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 显示图形
print(p2)




#########
#各模块内不同组物种丰度的差异比较（本示例图省事儿直接使用 Tukey test，实际情况中慎重选择统计检验方法！）
# 加载所需的库
library(dplyr)
library(multcomp)

# 选择特定模块的数据
module_data <- filter(otu, module == "1")

# 对该模块进行ANOVA
fit <- aov(module_z_score ~ group, data = module_data)

# 输出ANOVA的结果
summary(fit)

# 进行Tukey HSD测试
tukey <- TukeyHSD(fit, "group")

# 查看Tukey HSD测试的结果
tukey
# 加载所需的库
library(dplyr)
library(multcomp)
library(ggplot2)

# 假设你已经有了一个名为otu的数据框，以及进行ANOVA和Tukey HSD测试的结果

# 选择特定模块的数据
module_data <- filter(otu, module == "7")

# 对该模块进行ANOVA
fit <- aov(module_z_score ~ group, data = module_data)

# 进行Tukey HSD测试
tukey <- TukeyHSD(fit, "group")

# 使用cld获取带有字母标注的比较结果
tukey_cld <- cld(glht(fit, linfct = mcp(group = "Tukey")), alpha = 0.05)

# 查看结果
print(tukey_cld)


