# Readme



## 直接运行程序

编译命令：

```plain
g++ main.cpp -o main -std=c++20
g++ main2.cpp -o main2 -std=c++20
```

运行命令：

```plain
./main <methodname> <s> <n> <type> <filename>
```

<methodname> 为方法名称。

<s> 为阶数。

<n> 为步数。即设置步长为 $\dfrac Tn$。

<type> 为0或1，表示作课件中的第一条曲线还是第二条曲线。

<filename> 为输出文件名，将解的所有数据输出到该文件中。此项可以没有。

## 用 Makefile 运行程序

`make story` 生成报告；

`make run` 生成报告除最后一节外所有表格所用的数据。

