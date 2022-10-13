// look_md_results.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "filedlg.h"
#include "readfile.h"


void Statistika (vector<double> & dF, const char * name)
{
	int rows = dF.size();
	vector<double> dF_means; 
	int N = 25;
	int wind = int(rows / N);

	double mat_ozh = 0.0;
	
	for (int j = 1; j <= N; j++)
	{
		double m = 0.0;
		int j_wind = j*wind;
		for (int i = (j-1)*wind; i < j_wind; i++)
		{
			m += dF[i];
		}
		m /= wind;
		dF_means.push_back(m);
		mat_ozh += m;
	}
	mat_ozh /= N;
	double D = 0.0;

	for (j = 0; j < N; j++)
	{
		double d = dF_means[j] - mat_ozh;
		D += d*d;
	}
	D /= N - 1;
	
	double sigma = sqrt(D / N);

	printf("\nStatistika %s N = %d\n", name, N);
	printf("mat_ozh = %f\n", mat_ozh);
	printf("sigma = %f\n", sigma);
}

#define WRITE_LOCKED_FORCES 0
#define WRITE_WORKED_FORCES 1

int main(int argc, char* argv[])
{
	TCHAR filter[] =     TEXT("Ghemical MD results File (*.txt)\0*.txt\0")
						 TEXT("All Files (*.*)\0*.*\0");
	TCHAR fpath[1024];
	TCHAR filename[1024];

	sprintf(filename, "\0");
	{
		DWORD nFilterIndex;

		vector<string> names;
		vector<string> *pnames = &names;
		vector<vector<double> > vectors;
		vectors.reserve(2000000);

		while (OpenFileDlg(0, filter, fpath, nFilterIndex) == S_OK)
		{		
			ReadDatFile(NULL, fpath, filename, &vectors, pnames);
			pnames = NULL;

			printf("\nfilename %s\n\n", filename);

			int cols = names.size();
			int rows = vectors.size();
			
#if WRITE_LOCKED_FORCES
			int cMom = 4 - 1;
			int cVx = 5 - 1;
			int cFxup = 14 - 1;
			int cFxdw = 17 - 1;

			int cVxup = 8 - 1;
			int cVxdw = 11 - 1;
#endif
#if WRITE_WORKED_FORCES
			int cMom = 4 - 1;

			int cVx = 5 - 1;
			int cVxup = 14 - 1;
			int cVxdw = 17 - 1;

			int cVx_wk_up = 8 - 1;
			int cVx_wk_dw = 11 - 1;

			int cFx_wk_up = 20 - 1;
			int cFx_wk_dw = 23 - 1;

#endif

			vector<double> means(cols, 0.0);


			printf("vectors.size() = %d\n",rows);
			printf("names.size() = %d\n", cols);

			for (vector<vector<double> >::iterator it = vectors.begin();
			it != vectors.end(); it++)
			{
				for (int c = 0; c < cols; c++)
				{
					means[c] += (*it).operator [](c);
				}
			}

			for (int c = 0; c < cols; c++)
			{
				means[c] /= rows;
				printf("mean(%s) = %f\n", names[c].c_str(), means[c]);
			}

			int r0 = 0;

			cout << "enter r0\n";
			cin >> r0;

#if WRITE_LOCKED_FORCES
			vector<double> dF(rows-r0);
			for (int r = r0; r < rows; r++)
			{
				dF[r-r0] = vectors[r][cFxup] - vectors[r][cFxdw];
			}

			Statistika (dF, "dF");

			vector<double> Mom(rows-r0);
			for (r = r0; r < rows; r++)
			{
				Mom[r-r0] = vectors[r][cMom];
			}

			Statistika (Mom, "Mom");

			vector<double> dV(rows-r0);
			for (r = r0; r < rows; r++)
			{
				dV[r-r0] = vectors[r][cVxup] - vectors[r][cVxdw];
			}

			Statistika (dV, "dV");

			vector<double> Vx(rows-r0);
			for (r = r0; r < rows; r++)
			{
				Vx[r-r0] = vectors[r][cVx];
			}

			Statistika (Vx, "Vx");
#endif
#if WRITE_WORKED_FORCES
			vector<double> dF_wk(rows-r0);
			for (int r = r0; r < rows; r++)
			{
				dF_wk[r-r0] = vectors[r][cFx_wk_up] - vectors[r][cFx_wk_dw];
			}

			Statistika (dF_wk, "dF_wk");


			vector<double> dV_wk(rows-r0);
			for (r = r0; r < rows; r++)
			{
				dV_wk[r-r0] = vectors[r][cVx_wk_up] - vectors[r][cVx_wk_dw];
			}

			Statistika (dV_wk, "dV_wk");

			//if (!worked[n1])
			vector<double> Mom(rows-r0);
			for (r = r0; r < rows; r++)
			{
				Mom[r-r0] = vectors[r][cMom];
			}

			Statistika (Mom, "Mom");

			vector<double> dV(rows-r0);
			for (r = r0; r < rows; r++)
			{
				dV[r-r0] = vectors[r][cVxup] - vectors[r][cVxdw];
			}

			Statistika (dV, "dV");

			vector<double> Vx(rows-r0);
			for (r = r0; r < rows; r++)
			{
				Vx[r-r0] = vectors[r][cVx];
			}

			Statistika (Vx, "Vx");
#endif
		}
	}
	/*else
	{
		DWORD nFilterIndex;
		if (SaveFileDlg(0, filename, filter, nFilterIndex) == S_OK)
		{
			SetDlgItemText(ref->hDlg,IDC_EDIT_TRAJFILE2, filename);
		}	
	}*/
	
	printf("Hello World!\n");
	return 0;
}

