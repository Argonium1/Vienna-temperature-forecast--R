# 加载包
library(TSA)
library(csv)
library(stats)
library(astsa)
library(forecast)
library(tseries)
library(urca)

# 读取CSV文件
temperature_data <- read.csv("F:\\Desktop\\时间序列分析\\Temperature_Vienna_OnlyCelsius.csv")
# 合并 month, day, year 列
temperature_data$date <- as.Date(paste(temperature_data$Year,temperature_data$Month, temperature_data$Day, sep = "-"))
# 删除原始的 month, day, year 列
temperature_data <- subset(temperature_data, select = -c(Month, Day, Year))
# 显示更新后的数据集的前几行
head(temperature_data)

# 筛选出所有值为-99的行
filtered_data_1 <- temperature_data[temperature_data$Celsius_Temperature == -72.77777778, ]
# 显示筛选后的数据集
print(filtered_data_1)
# 26条缺失值，手动补数据
temperature_data[835, "Celsius_Temperature"] =-2.740740741
temperature_data[836, "Celsius_Temperature"] =-4.481481481
temperature_data[841, "Celsius_Temperature"]=-3.888888889
temperature_data[842, "Celsius_Temperature"]=-5.166666667
temperature_data[852, "Celsius_Temperature"]=3.138888889
temperature_data[2107, "Celsius_Temperature"]=23.1
temperature_data[2108, "Celsius_Temperature"]=23.4222222
temperature_data[2109, "Celsius_Temperature"]=23.7444444
temperature_data[2110, "Celsius_Temperature"]=24.0666667
temperature_data[4004, "Celsius_Temperature"]=19.1388889
temperature_data[4397, "Celsius_Temperature"]=11.66666667
temperature_data[4594, "Celsius_Temperature"]=15.83333333
temperature_data[6358, "Celsius_Temperature"]=1.083333333
temperature_data[7051, "Celsius_Temperature"]=-0.407407407
temperature_data[7052, "Celsius_Temperature"]=-1.314814815
temperature_data[7178, "Celsius_Temperature"]=12.5
temperature_data[7317, "Celsius_Temperature"]=13.72222222
temperature_data[7368, "Celsius_Temperature"]=2.444444444
temperature_data[7787, "Celsius_Temperature"]=5
temperature_data[8101, "Celsius_Temperature"]=7.166666667
temperature_data[8102, "Celsius_Temperature"]=5
temperature_data[8262, "Celsius_Temperature"]=14.97222222
temperature_data[8284, "Celsius_Temperature"]=8.32222222
temperature_data[8285, "Celsius_Temperature"]=9.97777778
temperature_data[8286, "Celsius_Temperature"]=11.63333333
temperature_data[8287, "Celsius_Temperature"]=13.28888889
filtered_data_2 <- temperature_data[temperature_data$Celsius_Temperature == -72.77777778, ]
# 显示筛选后的数据集
print(filtered_data_2)

###至此，数据无错误

# 分训练集，测试集
# 创建训练集（2019-01-01之前）和测试集（2019-01-01之后）
train_data <- subset(temperature_data, date < "2019-01-01")
test_data <- subset(temperature_data, date >= "2019-01-01")
# 显示训练集和测试集的行数
cat("训练集观测数量：", nrow(train_data), "\n")
cat("测试集观测数量：", nrow(test_data), "\n")

#画原始数据
plot(train_data$date,train_data$Celsius_Temperature ,pch =16 , col = "RED", cex = 0.2)

# 将数据集转换为时间序列对象
temperature_ts <- ts(train_data$Celsius_Temperature,start = c(1996,9,11),frequency = 365)
# 显示时间序列对象的摘要信息
print(summary(temperature_ts))

# 季节性分解
seasonal_decompose_m <- decompose(temperature_ts, type = c( "multiplicative"))
seasonal_decompose_a <- decompose(temperature_ts, type = c("additive"))
# 绘制季节分解图
pdf("seasonal.pdf",width = 32, height = 4)
plot(seasonal_decompose_m)
plot(seasonal_decompose_a)
#明显的周期性

# 进行ADF检验
adf_test <- adf.test(temperature_ts)
print(adf_test)
# p值远远小于0.05，表示拒绝原假设，即时间序列是平稳的

# 绘制自相关函数（ACF）图
pdf("acf_365.pdf", width = 6, height = 8)
acf(temperature_ts, lag = 365 * 20, main = "Autocorrelation Plot")
#明显的周期性.

# 绘制偏自相关函数（PACF）图
pdf("pacf_365.pdf", width = 32*4, height = 8)
#图要展开画，要不然就是一坨黑
pacf(temperature_ts, lag = 365 * 20, main = "Partial Autocorrelation Plot")
#明显的周期性.

#数据在ADF检验时已经确认平稳，试试一次差分什么情况
# 对温度时间序列数据进行一次差分
diff_temperature <- diff(temperature_ts)
# 绘制差分后的时间序列图
plot(diff_temperature, main = "First Difference of Temperature Time Series", ylab = "Differenced Temperature")
# 进行ADF检验
adf_test <- adf.test(diff_temperature)
print(adf_test) 
# 绘制自相关函数（ACF）图
acf(diff_temperature, lag = 365 * 20, main = "Autocorrelation Plot")
# 绘制偏自相关函数（PACF）图
pacf(diff_temperature, main = "Partial Autocorrelation Plot")
#明显周期性

#看看趋势如何
#线性趋势
# 执行线性趋势分析
trend_lm <- lm(temperature_ts ~ time(temperature_ts))
# 显示回归模型的摘要信息
summary(trend_lm)
# 绘制数据和趋势线
plot(temperature_ts, main = "Temperature Time Series with Linear Trend", ylab = "Temperature (Celsius)")
abline(trend_lm, col = "red")

#余弦趋势模拟？
library(TSA)
har=harmonic(temperature_ts,1)
ml=lm(temperature_ts~har)
summary(ml)
pdf("余弦趋势模拟.pdf", width = 32, height = 4)
plot (ts (fitted(ml), freq=365,start=c(1996,9,11)),
      ylab= 'Temperature', type='l',
      ylim=range (c(fitted (ml) , temperature_ts)),col="red" ) ; points (temperature_ts ,pch =16 , col = "blue", cex = 0.2)
#效果还行？
#残差-时间图
pdf("残差-时间图.pdf", width = 32, height = 4)
plot (y=rstudent (ml), x=as.vector (time (temperature_ts) ) ,
      xlab= 'Time', ylab='Standardized Residuals', type='o',pch =16 , col = "blue", cex = 0.2)
#±4°C的差距还是比较大
#残差是否符合正态
hist(rstudent(ml),xlab="标准残差")
qqnorm(rstudent(ml))
#拟合一下
plot(forecast(ml))

#ARMA

# 绘制自相关函数（ACF）图
pdf("acf_365.pdf", width = 6, height = 8)
acf(temperature_ts, lag = 365 * 20, main = "Autocorrelation Plot")
#明显的周期性.

# 绘制偏自相关函数（PACF）图
pdf("pacf_365.pdf", width = 32*4, height = 8)
#图要展开画，要不然就是一坨黑
pacf(temperature_ts, lag = 365 * 20, main = "Partial Autocorrelation Plot")
#明显的周期性.
#AR(3),diff(0),MA(0)

#扩展自相关
eacf(temperature_ts)
#看不出来啊

#最优子集(这个函数运行太慢了，我手动进行一下)
library(forecast)
best_aic <- Inf
best_order <- c(0, 0, 0)
# 迭代尝试各种参数组合
for (p in 0:3) {
  for (d in 0:2) {
    for (q in 0:2) {
      # 拟合ARIMA模型
      model <- Arima(temperature_ts, order = c(p, d, q), seasonal = list(order = c(0, 1, 0)))
      
      # 计算AIC
      aic <- AIC(model)
      
      # 更新最佳参数和AIC
      if (aic < best_aic) {
        best_aic <- aic
        best_order <- c(p, d, q)
      }
    }
  }
}
# 打印最佳参数
print(best_order)

#自动arima，aic最小认为合适
#library(forecast)
best_fit <- auto.arima(temperature_ts, ic = "aic", seasonal = TRUE)
summary(best_fit)

#手动操作，d=0已经确定，acf说ar应该在2,3附近
arima_model <- arima(temperature_ts, order = c(2,0,0))
summary(arima_model)
arima_model <- arima(temperature_ts, order = c(2,0,1))
summary(arima_model)
arima_model <- arima(temperature_ts, order = c(2,0,2))
summary(arima_model)
arima_model <- arima(temperature_ts, order = c(2,0,3))
summary(arima_model)
arima_model <- arima(temperature_ts, order = c(3,0,0))
summary(arima_model)
arima_model <- arima(temperature_ts, order = c(3,0,1))
summary(arima_model)
arima_model <- arima(temperature_ts, order = c(3,0,2))
summary(arima_model)
arima_model <- arima(temperature_ts, order = c(3,0,3))
summary(arima_model)
#看下来是（2,0,3）

#季节ARIMA
#一次差分
temperature_ts_diff <- diff(temperature_ts)

# 打印差分后的前几个观测值
print(head(temperature_ts_diff))

plot(ts_diff)
acf(ts_diff,lag=365*20)
pacf(ts_diff,lag=365*20)

arima_model <- arima(temperature_ts, order = c(2,1,1),seasonal=list(c(0,1,0)))

# 使用拟合的ARIMA模型进行未来499个值的预测
forecast_values <- forecast(arima_model, h = 499)

# 加载 test_data 数据集
# 请确保你已经加载了 test_data 数据集，以便进行比较

# 绘制预测值和测试数据的对比图
plot(forecast_values)

lines(test_data, col = "red")

# 计算残差
residuals <- residuals(arima_model)

# 绘制残差-时间图
plot(residuals, xlab = "Time", ylab = "Residuals", type = "o", pch = 16, col = "blue")

# 绘制QQ图
qqnorm(residuals)
qqline(residuals)