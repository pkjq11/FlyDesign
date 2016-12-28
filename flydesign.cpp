#include "flydesign.h"


FlyDesign::FlyDesign()
{
}

double FlyDesign::D2R(double deg){
return deg / 180.0 * M_PI;
}

Point3d<double> FlyDesign::WGS2Mercator(Point3d<double> inPoint){
    Point3d<double> result;
    result.z = inPoint.z;
    double const MaxSize = 20037508.900383711000000;
    result.x = inPoint.x * MaxSize / 180;
    double y = log(tan((90 + inPoint.y) * M_PI / 360)) / (M_PI / 180);
    result.y= y * MaxSize / 180;
    return result;
}

Point3d<double> FlyDesign::Mercator2WGS(Point3d<double> inPoint){
    Point3d<double> result;
    result.z = inPoint.z;
    double const MaxSize = 20037508.900383711000000;
    result.x = inPoint.x / MaxSize * 180;
    double y = inPoint.y / MaxSize * 180;
    result.y = 180 / M_PI * (2 * std::atan(std::exp(y * M_PI / 180)) - M_PI / 2);
    return result;
}

bool FlyDesign::Isfit()
{
    if((H>H_ky_high)||H<H_ky_low)
        return false;
    else return true;
}

bool FlyDesign::DesignLine(const int& csfs_,const double& fi_,const double& H_ky_low_,const double& H_ky_high_,const double& H_,
                           vector<vector<Point3d<double> > > &Air_line,vector<vector<Point3d<double> > > &Direct_line)
{
    //参数赋值
    csfs=csfs_;
    fi=fi_;
    H_ky_low=H_ky_low_;
    H_ky_high=H_ky_high_;
    H=H_;

    if(In_points.size()<1||!Isfit())
    {
        flag=false;
        return false;
    }

    L=H/tan(D2R(fi))*csfs;
    if(In_points.size()==1)
    {
        double theta=0;
        vector<Point3d<double> > air_temp;
        Point3d<double> in_mercator = WGS2Mercator(In_points[0]);
        while(theta<=2*M_PI)
        {
            Point3d<double> air_mercator;
            air_mercator.y=in_mercator.y+(H/tan(D2R(fi)))*sin(theta);
            air_mercator.x=in_mercator.x+(H/tan(D2R(fi)))*cos(theta);
            air_mercator.z=in_mercator.z+H;
            air_temp.push_back(Mercator2WGS(air_mercator));
            theta=theta+M_PI/50;
        }
        Air_line.push_back(air_temp);
        Direct_line.push_back(air_temp);
    }
    else{
        //此次更新只考虑每条航线的首尾两点进行作图
        for(size_t i=1;i<In_points.size();i++)
        {
            vector<Point3d<double> > air_temp;
            vector<Point3d<double> > direct_temp;
            direct_temp.push_back(WGS2Mercator(In_points[i-1]));
            direct_temp.push_back(WGS2Mercator(In_points[i]));
            double deta_Y,deta_X;
            if(direct_temp[i].y==direct_temp[i-1].y&&direct_temp[i].x!=direct_temp[i-1].x)
                deta_Y=L,deta_X=0;
            else if(direct_temp[i].x==direct_temp[i-1].x&&direct_temp[i].y!=direct_temp[i-1].y)
                deta_Y=0,deta_X = L;
            else if(direct_temp[i].x==direct_temp[i-1].x&&direct_temp[i].y==direct_temp[i-1].y)
            {
                flag=false;
                return false;
            }
            else
            {
                double tanAlpha = (direct_temp[i].y-direct_temp[i-1].y)/(direct_temp[i].x-direct_temp[i-1].x);
                deta_Y=abs( L*cos(atan(tanAlpha)))*csfs;
                deta_X=abs( L*sin(atan(tanAlpha)))*csfs;
            }
            int xd=direct_temp[i-1].x<direct_temp[i].x?1:-1;
            int yd=direct_temp[i-1].y<direct_temp[i].y?1:-1;
            for(vector<Point3d<double> >::iterator iter=direct_temp.begin();iter!=direct_temp.end();iter++)
            {
                Point3d<double> P_air_mercator;
                P_air_mercator.y=(*iter).y+deta_Y*xd;
                P_air_mercator.x=(*iter).x-deta_X*yd;
                P_air_mercator.z=(*iter).z;
                air_temp.push_back(Mercator2WGS(P_air_mercator));
            }
            Air_line.push_back(air_temp);
            Direct_line.push_back(direct_temp);
        }
    }
    flag=true;
    return true;
}
bool FlyDesign::OutCoverBeam(const double& sita_,vector<vector<Point3d<double> > > &Direct_line,vector<vector<Point3d<double> > > &Cover_beam_near,vector<vector<Point3d<double> > > &Cover_beam_far)
{
    if(!flag||In_points.size()<=1)
    {
        flag=false;
        return false;
    }
    sita=sita_;
    double L1 = (H*(tan(D2R(90-fi))-tan(D2R(90-fi-0.5*sita))));
    double L2 = (H*(tan(D2R(90-fi+0.5*sita))-tan(D2R(90-fi))));
    for(vector<vector<Point3d<double> > >::iterator iter=Direct_line.begin();iter!=Direct_line.end();iter++)
    {
        double deta_Y1,deta_Y2,deta_X1,deta_X2;
        vector<Point3d<double> >near_line(2);
        vector<Point3d<double> >far_line(2);
        int xd=(*iter)[0].x<(*iter)[1].x?1:-1;
        int yd=(*iter)[0].y<(*iter)[1].y?1:-1;
        if((*iter)[1].y==(*iter)[0].y&&(*iter)[1].x!=(*iter)[0].x)
        {
            deta_Y1=L1*csfs;
            deta_Y2=L2*csfs;
            deta_X1=deta_X2=0;
        }
        else if((*iter)[1].y!=(*iter)[0].y&&(*iter)[1].x==(*iter)[0].x)
        {
            deta_X1=L1*csfs;
            deta_X2=L2*csfs;
            deta_Y1=deta_Y2=0;
        }
        else if((*iter)[1].y!=(*iter)[0].y&&(*iter)[1].x!=(*iter)[0].x)
        {
            double tanAlpha = ((*iter)[1].y-(*iter)[0].y)/((*iter)[1].x-(*iter)[0].x);
            deta_Y1=abs(L1*cos(atan(tanAlpha)))*csfs;
            deta_Y2=abs(L2*cos(atan(tanAlpha)))*csfs;
            deta_X1=abs(L1*sin(atan(tanAlpha)))*csfs;
            deta_X2=abs(L2*sin(atan(tanAlpha)))*csfs;
        }
        else
        {
            flag = false;
            return false;
        }
        near_line[0].x = (*iter)[0].x-yd*deta_X1;
        near_line[0].y = (*iter)[0].y+xd*deta_Y1;
        near_line[0].z = (*iter)[0].z;
        near_line[1].x = (*iter)[1].x-yd*deta_X1;
        near_line[1].y = (*iter)[1].y+xd*deta_Y1;
        near_line[1].z = (*iter)[1].z;
        far_line[0].x = (*iter)[0].x+yd*deta_X2;
        far_line[0].y = (*iter)[0].y-xd*deta_Y2;
        far_line[0].z = (*iter)[0].z;
        far_line[1].x = (*iter)[1].x+yd*deta_X2;
        far_line[1].y = (*iter)[1].y-xd*deta_Y2;
        far_line[1].z = (*iter)[1].z;
        near_line[0] = Mercator2WGS(near_line[0]);
        near_line[1] = Mercator2WGS(near_line[1]);
        far_line[0] = Mercator2WGS(far_line[0]);
        far_line[1] = Mercator2WGS(far_line[1]);
        Cover_beam_near.push_back(near_line);
        Cover_beam_far.push_back(far_line);
    }
    flag=true;
    return true;
}

bool FlyDesign::OutSwath(const double& AD_start_,const double& AD_length_,const double& fs_,const double& tao_,
                         vector<vector<Point3d<double> > > &Direct_line,vector<vector<Point3d<double> > > &Swath_start,vector<vector<Point3d<double> > > &Swath_end)
{
    if(!flag||In_points.size()<=1)
    {
        flag=false;
        return false;
    }
    AD_start=AD_start_;
    AD_length=AD_length_;
    fs=fs_;
    tao=tao_;
    double W=((AD_length*1000/fs-tao)*300);
    deta_M=abs(H/tan(D2R(fi)))-sqrt(pow(AD_start*300,2)-H*H);
    if(deta_M<0||W<0||W<deta_M)
      {
        flag=false;
        return false;
      }
    double deta_N=deta_M-W;//波束覆盖另一侧的相反数
    for(vector<vector<Point3d<double> > >::iterator iter=Direct_line.begin();iter!=Direct_line.end();iter++)
    {
        double deta_Y1,deta_Y2,deta_X1,deta_X2;
        vector<Point3d<double> >start_line(2);
        vector<Point3d<double> >end_line(2);
        int xd=(*iter)[0].x<(*iter)[1].x?1:-1;
        int yd=(*iter)[0].y<(*iter)[1].y?1:-1;
        if((*iter)[1].y==(*iter)[0].y&&(*iter)[1].x!=(*iter)[0].x)
        {
            deta_Y1=deta_M*csfs;
            deta_Y2=deta_N*csfs;
            deta_X1=deta_X2=0;
        }
        else if((*iter)[1].y!=(*iter)[0].y&&(*iter)[1].x==(*iter)[0].x)
        {
            deta_X1=deta_M*csfs;
            deta_X2=deta_N*csfs;
            deta_Y1=deta_Y2=0;
        }
        else if((*iter)[1].y!=(*iter)[0].y&&(*iter)[1].x!=(*iter)[0].x)
        {
            double tanAlpha = ((*iter)[1].y-(*iter)[0].y)/((*iter)[1].x-(*iter)[0].x);
            deta_Y1=abs(deta_M*cos(atan(tanAlpha)))*csfs;
            deta_Y2=abs(deta_N*cos(atan(tanAlpha)))*csfs;
            deta_X1=abs(deta_M*sin(atan(tanAlpha)))*csfs;
            deta_X2=abs(deta_N*sin(atan(tanAlpha)))*csfs;
        }
        else
        {
            flag = false;
            return false;
        }
        start_line[0].x = (*iter)[0].x-yd*deta_X1;
        start_line[0].y = (*iter)[0].y+xd*deta_Y1;
        start_line[0].z = (*iter)[0].z;
        start_line[1].x = (*iter)[1].x-yd*deta_X1;
        start_line[1].y = (*iter)[1].y+xd*deta_Y1;
        start_line[1].z = (*iter)[1].z;
        end_line[0].x = (*iter)[0].x+yd*deta_X2;
        end_line[0].y = (*iter)[0].y-xd*deta_Y2;
        end_line[0].z = (*iter)[0].z;
        end_line[1].x = (*iter)[1].x+yd*deta_X2;
        end_line[1].y = (*iter)[1].y-xd*deta_Y2;
        end_line[1].z = (*iter)[1].z;
        start_line[0] = Mercator2WGS(start_line[0]);
        start_line[1] = Mercator2WGS(start_line[1]);
        end_line[0] = Mercator2WGS(end_line[0]);
        end_line[1] = Mercator2WGS(end_line[1]);
        Swath_start.push_back(start_line);
        Swath_end.push_back(end_line);

    }
    flag=true;
    return true;

}

bool FlyDesign::OutNESZ(const double& LightVs_,const double& Fre_,const double& gain_,const double& Fn_,const double& Ls_,
                        const double& Va_,const double& NESZ_YQ_,const double& Ka_,const double& Kr_,const double& k_,
                        const double& T0_,const double& Res_,const double& TpRatio_,const double& Pfzh_,
                        vector<P_NESZ> &R_NESZ,vector<P_NESZ> &R_NESZ_YQ)
{
    if(!flag)
    {
        flag=false;
        return false;
    }
    LightVs=LightVs_;
    Fre=Fre_;
    gain=gain_;
    Fn=Fn_;
    Ls=Ls_;
    Va=Va_;
    NESZ_YQ=NESZ_YQ_;
    Ka=Ka_;
    Kr=Kr_;
    k=k_;
    T0=T0_;
    Res=Res_;
    TpRatio=TpRatio_;
    Pfzh=Pfzh_;
    double lamda=LightVs/Fre;
    double GtGr=pow(10,(gain/10))*pow(10,(gain/10));
    double temp1=2*pow(4*M_PI,3)*k*T0*pow(10,Fn/10)*pow(10,Ls/10)*(Va/3.6);
    double temp2=GtGr*pow(lamda,3)*Ka*Kr*Res*Pfzh*TpRatio;
    double R_min=sqrt((L*csfs-deta_M)*(L*csfs-deta_M)+H*H);
    double R_max=sqrt((L*csfs+deta_M)*(L*csfs+deta_M)+H*H);
    double R_temp=R_min;
    while(R_temp<=R_max)
    {
       double NESZ=10*log10(temp1*pow(R_temp,3)*cos(asin(H/R_temp))/temp2);
       R_NESZ.push_back(make_pair(R_temp,NESZ));
       R_NESZ_YQ.push_back(make_pair(R_temp,NESZ_YQ));
       R_temp=R_temp+1;
    }

    return true;

}
bool FlyDesign::BackStep(vector<Point3d<double> > &Air_line,vector<Point3d<double> > &Original_line)
{
    Point3d<double> air_mercator1 = WGS2Mercator(Air_line[0]);
    Point3d<double> air_mercator2 = WGS2Mercator(Air_line[1]);
    Point3d<double> P_original1,P_original2;
    double deta_Y,deta_X;
    if(air_mercator2.y==air_mercator1.y&&air_mercator2.x!=air_mercator1.x)
        deta_Y=L*csfs,deta_X=0;
    else if(air_mercator2.x==air_mercator1.x&&air_mercator2.y!=air_mercator1.y)
        deta_Y=0,deta_X = L*csfs;
    else if(air_mercator2.x==air_mercator1.x&&air_mercator2.y==air_mercator1.y)
    {
        flag=false;
        return false;
    }
    else
    {
        double tanAlpha = (air_mercator2.y-air_mercator1.y)/(air_mercator2.x-air_mercator1.x);
        deta_Y=abs( L*cos(atan(tanAlpha)));
        deta_X=abs( L*sin(atan(tanAlpha)));
    }
    int xd=air_mercator1.x<air_mercator2.x?1:-1;
    int yd=air_mercator1.y<air_mercator2.y?1:-1;
    P_original1.x = air_mercator1.x-yd*deta_X;
    P_original1.y = air_mercator1.y+xd*deta_Y;
    P_original1.z = air_mercator1.z;
    P_original2.x = air_mercator2.x-yd*deta_X;
    P_original2.y = air_mercator2.y+xd*deta_Y;
    P_original2.z = air_mercator1.z;
    P_original1 = Mercator2WGS(P_original1);
    P_original2 = Mercator2WGS(P_original2);
    Original_line.push_back(P_original1);
    Original_line.push_back(P_original2);
    return true;
}
