#ifndef _SOLVE_H_
#define _SOLVE_H_

#include <ctime>

#include <petscdm.h>
#include <petscdmda.h>

#include "progress_bar.hpp"

using namespace std;

class TransLine
{
public:
    Vec                     U = nullptr;                    // 电压、电流、电荷量Vector
    Vec                     I = nullptr;
    Vec                     Q = nullptr;
    
    DM                      da = nullptr;                   // 分布数据管理器

    PetscReal               I0;                             // 下一时刻左端电流
    PetscReal               In;                             // 下一时刻右端电流

    PetscReal               C0 = 200e-12;                     // 分布参数
    PetscReal               G0 = 0;
    PetscReal               L0 = 1e-6;
    PetscReal               R0 = 1;

    PetscReal               x0 = 0;                         // 迭代参数
    PetscReal               xL = 1;
    PetscReal               dx = 0.01;

    PetscReal               maxTime = 1e-6;
    PetscReal               t0 = 0;
    PetscReal               dt = 0.1*1e-9;

    PetscReal               realTime = 0;                   // 迭代控制参数
    PetscInt                realStep = 0;

    bool                    isStop = false;                 // 标志位
    bool                    isMonitor = false;

    PetscInt                maxSteps = -1;                  // 计算得到参数
    PetscInt                Mx = -1;

public:
    TransLine();

    ~TransLine();

    void init();                                            // 初始化

    void initState();                                       // 初始化U I Q

    void update(PetscReal U0, PetscReal Un);                // 一次迭代  U0, U1分别是左端、右端边界电压

    void monitor(string vecFileName="../vec.dat");          // 控制U,I,Q写入文件

    void stepping();                                        // 迭代次数+1

    void printInfo();                                       // 打印transSolve的信息
};

#endif