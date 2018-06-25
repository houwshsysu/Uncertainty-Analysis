using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using ESRI.ArcGIS;


namespace VECTOR_MODEL
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            if (!RuntimeManager.Bind(ProductCode.EngineOrDesktop))
            {
                MessageBox.Show("Unable to bind to ArcGIS runtime. Application will be shut down.");
                return;
            }

            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new A_Form1());

            #region CreateNodes
            //CreateNodes
            CreateNodes e = new CreateNodes();
            e.CreateNode();
            e.GetEveryKindSingleProbability();
            e.GetSurfaceNode();
            e.CreatFeatureClasses();
            e.CreatTINs();
            #endregion

            #region GetOriginal
            //GetOriginal
            GetOriginal g = new GetOriginal();
            g.GetTinData();
            g.GetSurfaceNormalVector();
            g.GetNodeNormalVector();
            g.GetNodeCurvature();
            g.CreateFeatureClass();
            g.CreateFeature();
            #endregion

            #region VectorModel
            //VectorModel
            VectorModel v = new VectorModel();
            for (int d = 0; d < VectorModel.TinCount; d++)
            {
                v.D = d;
                v.GetTinData(d);
                if (d == 0)
                {
                    //只需要一个
                    v.GetMeanInfo();
                }
                v.GetSurfaceNormalVector();
                v.GetNodeNormalVector();
                v.GetNodeCurvature();
                if (d == 0)
                {
                    //只需要一个
                    v.CreateFeatureClass();
                }
                v.CreateFeature();
            }
            #endregion

            //#region GravityVector
            ////GravityVector
            //GravityVector gv = new GravityVector();
            //for (int d = 0; d < VectorModel.TinCount; d++)
            //{
            //    gv.GetTinData(d);
            //    gv.GetVolume();
            //    gv.GetBarycentre();
            //    gv.CreatEnclose();
            //    if (d == 0)
            //    {
            //        //只需要一个
            //        gv.CreateFeatureClass();
            //    }
            //}
            //gv.CreateFeature();
            //#endregion
        }
    }
}

