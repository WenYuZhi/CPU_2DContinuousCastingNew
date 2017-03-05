#define Section 12
#define CoolSection 8
#define MoldSection 4
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void Physicial_Parameters(float);
float Boundary_Condition(int j, int ny, float Ly, float *, float *);
void Monitor_MeanTemperature(int, int, int *);
void OnestepSolve(int nx, int ny, float tao, float dx, float dy, float* H_Init);

float ***T;
float **T_Last;
float **T_New;
float **T_HoldLast;
float *Mean_TSurface, *Mean_TCentral;
float Vcast = -0.02, h = 1000.0, lamd = 50.0, Ce = 540.0, pho = 7000.0, a = 1.0;
float ccml[Section + 1] = { 0.0,0.2,0.4,0.6,0.8,1.0925,2.27,4.29,5.831,9.6065,13.6090,19.87014,28.599 };
float Taim[Section] = { 966.149841, 925.864746, 952.322083, 932.175537, 914.607117, 890.494263, 870.804443, 890.595825 };
float Lx = 0.125, Ly = 28.599, t_final = 1500.0, T_Cast = 1558, Tw = 30;

int main()
{
	FILE *fp = NULL;
	float **T_Predictive_Last;
	float **T_Predictive_New;
	int nx = 11, ny = 3001, tnpts = 7501, inter_num = 500, Num = 5;
	float tao, dx, dy, dh = 1.0, arf1, arf2, step = -0.002;
	float H_Init[Section] = { 1380,1170,980,800,1223.16,735.05,424.32,392.83,328.94,281.64,246.16,160.96 };
	float H_Init_Temp[Section] = { 1380,1170,980,800,1223.16,735.05,424.32,392.83,328.94,281.64,246.16,160.96 };
	float *y, **Mean_TSurfaceElement, **Mean_TSurfaceElementOne, *Delta_H_Init, **JacobianMatrix;
	bool disout = true;
	int Y_Label[Section + 1] = { 0 };

	clock_t begin, end, duration;

	y = (float*)calloc(ny, sizeof(float));
	Mean_TSurface = (float*)calloc(Section, sizeof(float));
	Mean_TCentral = (float*)calloc(Section, sizeof(float));
	Delta_H_Init = (float*)calloc(CoolSection, sizeof(float));

	T_New = (float**)calloc(nx, sizeof(float));
	for (int i = 0; i < nx; i++)
		T_New[i] = (float*)calloc(ny, sizeof(float));

	T_Last = (float**)calloc(nx, sizeof(float));
	for (int i = 0; i < nx; i++)
		T_Last[i] = (float*)calloc(ny, sizeof(float));

	T_HoldLast = (float**)calloc(nx, sizeof(float));
	for (int i = 0; i < nx; i++)
		T_HoldLast[i] = (float*)calloc(ny, sizeof(float));

	T_Predictive_New = (float**)calloc(nx, sizeof(float));
	for (int i = 0; i < nx; i++)
		T_Predictive_New[i] = (float*)calloc(ny, sizeof(float));

	T_Predictive_Last = (float**)calloc(nx, sizeof(float));
	for (int i = 0; i < nx; i++)
		T_Predictive_Last[i] = (float*)calloc(ny, sizeof(float));

	T = (float***)calloc(nx, sizeof(float));
	for (int i = 0; i < nx; i++)
		T[i] = (float**)calloc(ny, sizeof(float));
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			T[i][j] = (float*)calloc(21, sizeof(float));

	Mean_TSurfaceElement = (float**)calloc(CoolSection, sizeof(float));
	for (int i = 0; i < CoolSection; i++)
		Mean_TSurfaceElement[i] = (float*)calloc(CoolSection, sizeof(float));

	Mean_TSurfaceElementOne = (float**)calloc(CoolSection, sizeof(float));
	for (int i = 0; i < CoolSection; i++)
		Mean_TSurfaceElementOne[i] = (float*)calloc(CoolSection, sizeof(float));

	JacobianMatrix = (float**)calloc(CoolSection, sizeof(float));
	for (int i = 0; i < CoolSection; i++)
		JacobianMatrix[i] = (float*)calloc(CoolSection, sizeof(float));

	dx = Lx / (nx - 1);
	dy = Ly / (ny - 1);
	tao = t_final / (tnpts - 1);

	for (int j = 0; j < ny; j++)
		y[j] = j * dy;

	Y_Label[Section] = ny - 1;
	for (int j = 0; j < ny - 1; j++)
		for (int i = 1; i < Section; i++)
			if (y[j] <= ccml[i] && y[j + 1] >= ccml[i])
				Y_Label[i] = j + 1;

	printf("Casting Temperature = %f ", T_Cast);
	printf("\n");
	printf("The thick of steel billets(m) = %f ", Lx);
	printf("\n");
	printf("The length of steel billets(m) = %f ", Ly);
	printf("\n");
	printf("dx(m) = %f ", dx);
	printf("dy(m) = %f ", dy);
	printf("tao(s) = %f ", tao);
	printf("\n");
	printf("simulation time(s) = %f ", t_final);

	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			T_Last[i][j] = T_Cast;

	begin = clock();
	fp = fopen("C:\\Temperature2D_Static.txt", "a");
	for (int k = 0; k < tnpts - 1; k++)
	{
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				T_HoldLast[i][j] = T_Last[i][j];

		for (int m = 0; m < CoolSection + 1; m++)
		{
			if (m == CoolSection)
			{
				for (int PNum = 0; PNum < Num; PNum++)
				{
					for (int temp = 0; temp < Section; temp++)
						H_Init_Temp[temp] = H_Init[temp];
					OnestepSolve(nx, ny, tao, dx, dy, H_Init_Temp);
					for (int j = 0; j < ny; j++)
						for (int i = 0; i < nx; i++)
							T_Last[i][j] = T_New[i][j];
				} 

				Monitor_MeanTemperature(ny, nx, Y_Label);
				for (int temp = 0; temp < CoolSection; temp++)
				    for (int column = 0; column < CoolSection; column++)
					    Mean_TSurfaceElementOne[temp][column] = Mean_TSurface[column + MoldSection];
			}

			else
			{
				for (int PNum = 0; PNum < Num; PNum++) 
				{
					for (int temp = 0; temp < Section; temp++)
						H_Init_Temp[temp] = H_Init[temp];
					if( PNum == 0 )
					    H_Init_Temp[m + MoldSection] = H_Init[m + MoldSection] + dh;
					OnestepSolve(nx, ny, tao, dx, dy, H_Init_Temp);
					for (int j = 0; j < ny; j++)
						for (int i = 0; i < nx; i++)
							T_Last[i][j] = T_New[i][j];
				}
				
				Monitor_MeanTemperature(ny, nx, Y_Label);
				for (int column = 0; column < CoolSection; column++)
					Mean_TSurfaceElement[m][column] = Mean_TSurface[column + MoldSection];
			}

			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
					T_Last[i][j] = T_HoldLast[i][j];
		}

		if (k % Num == 0)
		{
			printf("\nJacobianMatrix=\n");
			for (int row = 0; row < CoolSection; row++)
			{
				for (int column = 0; column < CoolSection; column++)
				{
					JacobianMatrix[row][column] = (Mean_TSurfaceElement[row][column] - Mean_TSurfaceElementOne[row][column]) / dh;
					printf("%f ", JacobianMatrix[row][column]);
				}
				printf("\n");
			}

			for (int temp = 0; temp < CoolSection; temp++)
				Delta_H_Init[temp] = 0.0;

			printf("\nDelta_H_Init=\n");
			for (int temp = 0; temp < CoolSection; temp++)
			{
				for (int column = 0; column < CoolSection; column++)
					Delta_H_Init[temp] += (Mean_TSurfaceElementOne[temp][column] - Taim[column]) * JacobianMatrix[temp][column];
				printf(" %f, ", Delta_H_Init[temp]);
			}

			arf1 = 0.0, arf2 = 0.0;
			for (int temp = 0; temp < CoolSection; temp++)
			{
				for (int lab = 0; lab < CoolSection; lab++)
				{
					arf1 += (Mean_TSurfaceElementOne[0][temp] - Taim[temp]) * JacobianMatrix[temp][lab] * Delta_H_Init[lab];
					arf2 += JacobianMatrix[temp][lab] * Delta_H_Init[lab] * JacobianMatrix[temp][lab] * Delta_H_Init[lab];
				}
			}

				printf("\n");
				printf(" Output time step = %d", k);
				printf(" Simulation time = %f", k * tao);
				printf("\n step = %f, arf1 = %f, arf2 = %f\n", step, arf1, arf2);
				printf("\n");

			for (int temp = 0; temp < CoolSection; temp++)
				H_Init[temp + MoldSection] += step *(Delta_H_Init[temp]);
		}

		for (int temp = 0; temp < Section; temp++)
			H_Init_Temp[temp] = H_Init[temp];
		OnestepSolve(nx, ny, tao, dx, dy, H_Init_Temp);
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				T_Last[i][j] = T_New[i][j];

		if(k % Num == 0)
		{
			printf("\nTSurface=\n");
			for (int temp = 0; temp < CoolSection; temp++)
				printf("%f, ", Mean_TSurface[temp + MoldSection]);

			printf("\nTSurface - Taim=\n");
			for (int temp = 0; temp < CoolSection; temp++)
				printf("%f, ", Mean_TSurface[temp + MoldSection] - Taim[temp]);

			for (int column = MoldSection; column < Section; column++)
			{
				fprintf(fp, "%f ", Mean_TSurface[column]);
				fprintf(fp, "%f, ", H_Init[column]);
			}
			fprintf(fp, "\n");
		}
		
	}

	fclose(fp);
	end = clock();
	printf("\n");
	printf("running time = %d (mseconds)\n", (end - begin));

	return 0;
}


void Physicial_Parameters(float T)
{
	float Ts = 1462.0, Tl = 1518.0, lamds = 30, lamdl = 50, phos = 7000, phol = 7500, ce = 540.0, L = 265600.0, fs = 0.0;
	if (T < Ts)
	{
		fs = 0;
		pho = phos;
		lamd = lamds;
		Ce = ce;
	}

	if (T >= Ts&&T <= Tl)
	{
		fs = (T - Ts) / (Tl - Ts);
		pho = fs*phos + (1 - fs)*phol;
		lamd = fs*lamds + (1 - fs)*lamdl;
		Ce = ce + L / (Tl - Ts);
	}

	if (T > Tl)
	{
		fs = 1;
		pho = phol;
		lamd = lamdl;
		Ce = ce;
	}

}

float Boundary_Condition(int j, int ny, float Ly, float *ccml_zone, float *H_Init)
{
	float YLabel, h;
	YLabel = (j * Ly) / float(ny - 1);

	for (int i = 0; i < Section; i++)
	{
		if (YLabel >= *(ccml_zone + i) && YLabel <= *(ccml_zone + i + 1))
		{
			h = *(H_Init + i);
		}
	}
	return h;
}

void Monitor_MeanTemperature(int ny, int nx, int *Y_Label)
{

	float temp_surface[Section] = { 0.0 }, temp_central[Section] = { 0.0 };
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < Section; i++)
			if (j > *(Y_Label + i) && j <= *(Y_Label + i + 1))
			{
				temp_surface[i] = temp_surface[i] + T_New[0][j];
				temp_central[i] = temp_central[i] + T_New[nx - 1][j];
			}

	//printf("\n Surface temperature\n");

	for (int i = 0; i < Section; i++)
	{
		Mean_TSurface[i] = temp_surface[i] / (*(Y_Label + i + 1) - *(Y_Label + i));
		//printf("zone %d = %f  ", i + 1, Mean_TSurface[i]);
	}

	//printf("\n Central temperature\n");
	for (int i = 0; i < Section; i++)
	{
		Mean_TCentral[i] = temp_central[i] / (*(Y_Label + i + 1) - *(Y_Label + i));
		//printf("zone %d = %f  ", i + 1, Mean_TCentral[i]);
	}
}



void OnestepSolve(int nx, int ny, float tao, float dx, float dy, float* H_Init)
{
	float T_Up, T_Down, T_Right, T_Left;

		for (int j = 0; j < ny; j++)
		{
			h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
			for (int i = 0; i < nx; i++)
			{
				Physicial_Parameters(T_Last[i][j]);
				a = lamd / (pho*Ce);
				if (i == 0 && j != 0 && j != ny - 1)  //1
				{
					T_Up = T_Last[i + 1][j];
					T_Down = T_Last[i + 1][j] - 2 * dx*h*(T_Last[i][j] - Tw) / lamd;
					T_Right = T_Last[i][j + 1];
					T_Left = T_Last[i][j - 1];
					T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
						a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
				}
				else if (i == nx - 1 && j != 0 && j != ny - 1)//2
				{
					T_Up = T_Last[i - 1][j];
					T_Down = T_Last[i - 1][j];
					T_Right = T_Last[i][j + 1];
					T_Left = T_Last[i][j - 1];
					T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
						a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
				}
				else if (j == 0 && i != 0 && i != nx - 1)//3
				{
					T_New[i][j] = T_Cast;
				}
				else if (j == ny - 1 && i != 0 && i != nx - 1)//4
				{
					T_Up = T_Last[i + 1][j];
					T_Down = T_Last[i - 1][j];
					T_Right = T_Last[i][j - 1];
					T_Left = T_Last[i][j - 1];
					T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
						a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
				}
				else if (i == 0 && j == 0)//5
				{
					T_New[i][j] = T_Cast;
				}
				else if (i == 0 && j == ny - 1)//6
				{
					T_Up = T_Last[i + 1][j];
					T_Down = T_Last[i + 1][j] - 2 * dx*h*(T_Last[i][j] - Tw) / lamd;
					T_Right = T_Last[i + 1][j];
					T_Left = T_Last[i][j - 1];
					T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
						a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
				}
				else if (i == nx - 1 && j == 0)//7
				{
					T_New[i][j] = T_Cast;
				}
				else if (i == nx - 1 && j == ny - 1)//8
				{
					T_Up = T_Last[i - 1][j];
					T_Down = T_Last[i - 1][j];
					T_Right = T_Last[i - 1][j];
					T_Left = T_Last[i][j - 1];
					T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
						a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
				}
				else//9
				{
					T_Up = T_Last[i + 1][j];
					T_Down = T_Last[i - 1][j];
					T_Right = T_Last[i][j + 1];
					T_Left = T_Last[i][j - 1];
					T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
						a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
				}
			}
		}
	
		/*for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				T_Last[i][j] = T_New[i][j];*/
}