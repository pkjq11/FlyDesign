#ifndef FLYDESIGN_H
#define FLYDESIGN_H

#include "flydesign_global.h"
#include <vector>
#include <math.h>
#include <QDebug>
#define M_PI 3.14159265358979323846
using namespace std;
template <class Type> struct Point3d
{
    Type x;                                //经度，单位：度
    Type y;                                //纬度，单位：度
    Type z;                                //拔高度，单位：米
    Point3d(){x=0;y=0;z=0;}
};

typedef pair<double,double> P_NESZ ;

class FLYDESIGNSHARED_EXPORT FlyDesign
{

public:
    FlyDesign();

    bool Isfit();                     //判断参数是否合理，合理返回true，否则返回false
    bool DesignLine(const int& csfs_,const double& fi_,const double& H_ky_low_,const double& H_ky_high,const double& H_,
                    vector<vector<Point3d<double> > > &Air_line,vector<vector<Point3d<double> > > &Direct_line);                //规划航线
    bool OutCoverBeam(const double& sita_,vector<vector<Point3d<double> > > &Direct_line,vector<vector<Point3d<double> > > &Cover_beam_near,vector<vector<Point3d<double> > > &Cover_beam_far);              //计算波束覆盖区域
    bool OutSwath(const double& AD_start_,const double& AD_length_,const double& fs_,const double& tao_,
                  vector<vector<Point3d<double> > > &Direct_line,vector<vector<Point3d<double> > > &Swath_start,vector<vector<Point3d<double> > > &Swath_end);                  //计算测绘带
    bool OutNESZ(const double& LightVs_,const double& Fre_,const double& gain_,const double& Fn_,const double& Ls_,
                 const double& Va_,const double& NESZ_YQ_,const double& Ka_,const double& Kr_,const double& k_,
                 const double& T0_,const double& Res_,const double& TpRatio_,const double& Pfzh_,
                 vector<P_NESZ> &R_NESZ,vector<P_NESZ> &R_NESZ_YQ);                   //判断测绘带内NESZ值

public:
    vector<Point3d<double> > In_points;         //输入的地面目标点集
    Point3d<double> Mercator2WGS(Point3d<double> inPoint); //WGS84坐标系转为web墨卡托
    Point3d<double> WGS2Mercator(Point3d<double> inPoint); //web墨卡托转为WGS84

private:
    double D2R(double deg);            //角度值转弧度值
//***************************需要输入的参数********************************************
private:
    int csfs;                          //雷达侧视方式：1为右侧视，-1为左侧视；
    double fi;                         //雷达下视角，单位：度
    double H_ky_low;                   //载机准飞空域高度（下限），单位：米
    double H_ky_high;                  //载机准飞空域高度（上限），单位：米
    double H;                          //载机飞行高度，单位：米
    double sita;                       //雷达俯仰向波束宽度
    double AD_start;                   //采样起始，单位：us
    double AD_length;                  //采样深度，单位：Kbit
    double fs;                         //采样频率，单位：MHz
    double tao;                        //脉冲宽度，单位：us
    double LightVs;                    //光速，单位：m/s
    double Fre;                        //中心频率，单位：Hz
    double gain;                       //系统增益，单位：dB
    double Fn;                         //接收机噪声系数，单位：dB
    double Ls;                         //双程气象损耗+系统损耗，单位：dB
    double Va;                         //载机飞行速度，单位：km/h
    double NESZ_YQ;                    //NESZ要求
    double Ka;                         //ka为方位、距离脉压加窗损失（在海明窗加权时，选-1.34dB）
    double Kr;                         //kr为方位、距离脉压加窗损失（在海明窗加权时，选-1.34dB）
    double k;                          //玻尔兹曼常数
    double T0;                         //接收机绝对温度，以290度计
    double Res;                        //地距分辨率
    double TpRatio;                    //占空比
    double Pfzh;                       //峰值功率
//****************************中间计算结果***********************************************
private:
    double L;                          //航线与目标方向性直线间的地面投影距离
    double deta_M;
    bool flag;
};

#endif // FLYDESIGN_H
