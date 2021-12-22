#include "solve.hpp"


TransLine::TransLine()
{

}

TransLine::~TransLine()
{
    if (nullptr != this->U) 
    {
        VecDestroy(&this->U);
        this->U = nullptr;
    }
    if (nullptr != this->I)
    {
        VecDestroy(&this->I);
        this->I = nullptr;
    }
    if (nullptr != this->Q) 
    {
        VecDestroy(&this->Q);
        this->Q = nullptr;
    }
    if (nullptr != this->da) 
    {
        DMDestroy(&this->da);
        this->da = nullptr;
    }
}

void TransLine::init()
{
    // Ts参数设置
    this->maxSteps = (this->maxTime - this->t0) / this->dt + 0.5;
    this->Mx = (this->xL - this->x0) / this->dx + 1;

    // 迭代步数从0开始，realTime = realStep * dt
    this->realStep = 0;
    this->realTime = this->t0;
    this->isStop = false;

    // 构造
    if (nullptr == da)
    {
        DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, this->Mx, 1, 1, NULL, &this->da);
        DMSetFromOptions(this->da);
        DMSetUp(this->da);
    }

    // Vec
    if (nullptr == this->U && nullptr == this->I && nullptr == this->Q)
    {
        DMCreateGlobalVector(this->da, &this->U);
        VecDuplicate(this->U, &this->I);
        VecDuplicate(this->U, &this->Q);
    }

    this->printInfo();

    this->initState();
}

void TransLine::initState()
{
    // 设置x的初值
    VecSet(this->U, 0);
    VecSet(this->I, 0);
    VecSet(this->Q, 0);
}

void TransLine::update(PetscReal U0, PetscReal Un)
{
    PetscScalar             *uptr, *iptr, *qptr;                    // 数据指针
    PetscScalar             *unext, *inext, *qnext;
    Vec                     temU, temI, temQ;
    PetscReal               a11, a12, a13, a21, a22, a23, a31, a32; // 迭代系数
    PetscReal               b11, b12, b21, b22, b31;
    PetscReal               upm, unm, ipm, inm, qpm, qnm;

    // 系数
    a11 = 0.5;      a12 = -0.5*this->dt/this->dx/this->C0;      a13 = -0.25*this->dt*this->G0/this->C0;
    a21 = 0.5;      a22 = -0.5*this->dt/this->dx/this->L0;      a23 = -0.25*this->dt*this->R0/this->L0;
    a31 = 0.5;      a32 = 0.25*this->dt;

    b11 = -1*this->dt/this->dx/this->C0;        b12 = -0.5*this->dt*this->G0/this->C0;
    b21 = -1*this->dt/this->dx/this->L0;        b22 = -0.5*this->dt*this->R0/this->L0;
    b31 = 0.5*this->dt;

    // 创建本地Vec，可以跨进程访问数据
    DMCreateGlobalVector(this->da, &temU);
    DMCreateGlobalVector(this->da, &temI);
    DMCreateGlobalVector(this->da, &temQ);

    VecCopy(this->U, temU);
    VecCopy(this->I, temI);
    VecCopy(this->Q, temQ);

    // 获取Vec数据指针
    DMDAVecGetArray(this->da, temU, &uptr);
    DMDAVecGetArray(this->da, temI, &iptr);
    DMDAVecGetArray(this->da, temQ, &qptr);

    DMDAVecGetArray(this->da, this->U, &unext);
    DMDAVecGetArray(this->da, this->I, &inext);
    DMDAVecGetArray(this->da, this->Q, &qnext);

    // loop
    for (size_t i = 0; i < this->Mx; i++)
    {
        if (0 == i)                             // 始端边界
        {
            PetscReal u1, i0;

            // update U k
            unext[i] = U0;
            uptr[i] = U0;

            u1 = uptr[1];
            i0 = iptr[0];

            // update I k+1
            inext[i] = this->dt / this->dx / this->L0 * 
                        (U0 - u1 - (this->dx*this->R0 - this->dx*this->L0/this->dt)*i0);
            
            this->I0 = inext[i];
        }
        else if (this->Mx-1 == i)               // 终端边界
        {
            // update Un
            unext[i] = Un;
            uptr[i] = Un;

            // update U
            inext[i] = this->dt / this->dx / this->L0 * 
                        (uptr[i-1] - uptr[i] - (this->dx*this->R0 - this->dx*this->L0/this->dt)*iptr[i]);

            // update I
            this->In = inext[i];
        }
        else                                    // 内部
        {
            // k+0.5时刻数据
            upm = (a11 + a13) * (uptr[i+1] + uptr[i]) + a12 * (iptr[i+1] - iptr[i]);
            unm = (a11 + a13) * (uptr[i] + uptr[i-1]) + a12 * (iptr[i] - iptr[i-1]);

            ipm = (a21 + a23) * (iptr[i+1] + iptr[i]) + a22 * (uptr[i+1] - uptr[i]);
            inm = (a21 + a23) * (iptr[i] + iptr[i-1]) + a22 * (uptr[i] - uptr[i-1]);

            qpm = a31 * (qptr[i+1] + qptr[i]) + a32 * (iptr[i+1] + iptr[i]);
            qnm = a31 * (qptr[i] + qptr[i-1]) + a32 * (iptr[i] + iptr[i-1]);

            // k+1
            unext[i] = uptr[i] + b11 * (ipm - inm) + b12 * (upm + unm);

            inext[i] = iptr[i] + b21 * (upm - unm) + b22 * (ipm + inm);

            qnext[i] = qptr[i] + b31 * (ipm + inm);
        }
    }
    
    // 数据反存到Vec
    DMDAVecRestoreArray(this->da, temU, &uptr);
    DMDAVecRestoreArray(this->da, temI, &iptr);
    DMDAVecRestoreArray(this->da, temQ, &qptr);

    DMDAVecRestoreArray(this->da, this->U, &unext);
    DMDAVecRestoreArray(this->da, this->I, &inext);
    DMDAVecRestoreArray(this->da, this->Q, &qnext);

    VecDestroy(&temU);
    VecDestroy(&temI);
    VecDestroy(&temQ);
}

void TransLine::monitor(string vecFileName)
{
    if (this->isMonitor)
    {
        PetscInt                N = this->maxSteps;
        PetscInt                N_save = N / 400;           // 总共保存文件数
        PetscInt                N_cout = N / 50;            // 进度条

        if (this->realStep % N_save == 0 || this->realStep == N)
        {
            PetscViewer myViewer;
            PetscViewerASCIIOpen(PETSC_COMM_WORLD, vecFileName.c_str(), &myViewer);
            VecView(this->U, myViewer);
            VecView(this->I, myViewer);
            VecView(this->Q, myViewer);
        }

        progressBar(this->realStep, this->maxSteps);
    }
}

void TransLine::stepping()
{
    // 步进
    if (this->realStep < this->maxSteps)
    {
        this->realStep++;
        this->realTime += this->dt;
    }
    else
    {
        this->isStop = true;
    }
}

void TransLine::printInfo()
{
    fmt::print("/n----------TransLine parameters----------/n");
    fmt::print("R0:                {}\n", this->R0);
    fmt::print("G0:                {}\n", this->G0);
    fmt::print("C0:                {}\n", this->C0);
    fmt::print("L0:                {}\n", this->L0);
    fmt::print("Time step:         {}\n", this->dt);
    fmt::print("Time range:        {} to {}\n", this->t0, this->maxTime);
    fmt::print("TransLine step:    {}\n", this->dx);
    fmt::print("TransLine range:   {} to {} \n\n", this->x0, this->xL);
}
