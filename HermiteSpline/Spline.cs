using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace HermiteSpline
{
    internal class Spline
    {
        internal readonly int funcnum, linknum, param_count, link_count;
        internal readonly double[,] A;
        internal readonly double[] X;

        struct herm_params
        {
            internal int type;
            internal int param_count; //к-сть параметрів у одній ланці
            internal int link_count;  //кількість ланок
            internal IntPtr A;       //параметри сплайна (link_conut*param_count елементів)
            internal IntPtr X;		 //точки наближення (link_count+1 елементів)
            internal IntPtr A128;
            internal IntPtr X128;
        };
        herm_params hp;
        /*
        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern double shit(double[] shits, Int64 count);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern Int64 alloc();

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern Int64 fre(Int64 ptr);
        */
        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern herm_params _HermGen(Byte funcnum, Byte linknum, double a, double b, double nu);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern void _free(herm_params hp);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern double _HermiteSpline(herm_params hp, double x, Byte derivative = 0);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern double _Func(Byte funcnum, double x, Byte derivative = 0);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern double _MaxError(herm_params hp, Byte funcnum, double from, double to);



        public Spline(int funcnum, int linknum, double a, double b, double nu)
        {
            try
            {
                this.hp = _HermGen((Byte)funcnum, (Byte)linknum, a, b, nu);
            }
            catch(Exception ex)
            {
                System.Windows.Forms.MessageBox.Show(ex.Message);
                return;
            }
            if (hp.X == IntPtr.Zero || hp.A == IntPtr.Zero)
            {
                throw new Exception("Couldn't create spline");
            }
            this.funcnum = funcnum;
            this.linknum = linknum;
            this.link_count = this.hp.link_count;
            this.param_count = this.hp.param_count;
            this.A = new double[this.link_count, this.param_count];
            this.X = new double[this.link_count + 1];

            Marshal.Copy(hp.X, X, 0, link_count + 1);
            double[] tmp = new double[link_count * param_count];
            Marshal.Copy((IntPtr)hp.A, tmp, 0, link_count * param_count);
            for (int i = 0; i < link_count; i++)
            {
                for (int j = 0; j < param_count; j++)
                {
                    A[i, j] = tmp[i * param_count + j];
                }
            }
            
        }

        ~Spline()
        {
            _free(this.hp);
        }

        public double Eval(double x)
        {
            return _HermiteSpline(hp, x);
        }
        public double OriginFunc(double x)
        {
            return _Func((Byte)funcnum, x);
        }
        public double EvalDer(double x)
        {
            return _HermiteSpline(hp, x, 1);
        }
        public double OriginDer(double x)
        {
            return _Func((Byte)funcnum, x, 1);
        }
        public double MaxError(double from, double to)
        {
            return _MaxError(hp, (Byte)funcnum, from, to);
        }

    }
}
