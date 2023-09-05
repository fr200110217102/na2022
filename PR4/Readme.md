# Readme

## 一键运行

生成报告文档：

```plain
make report
```

生成所有热方程的数据：

```plain
make runheat
```

生成所有对流方程的数据：

```plain
make runadv
```

## 自定义运行

编译

```plain
g++ main1.cpp -o main1
g++ main2.cpp -o main2
```

热方程的求解程序

```plain
./main1 <method> <n> <m> <filename>
```

对流方程的求解程序

```plain
./main2 <method> <n> <m> <filename>
```

这里`method`是方法名称，$n$和$m$分别是时间步数和空间步数，`filename`是输出文件名。