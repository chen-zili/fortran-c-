#ifndef TEMP_H
#define TEMP_H

#ifdef __cplusplus    //__cplusplus是cpp中自定义的一个宏
extern "C" {          //告诉编译器，这部分代码按C语言的格式进行编译，而不是C++的
#endif

    int temp(int n);

#ifdef __cplusplus
}
#endif

#endif