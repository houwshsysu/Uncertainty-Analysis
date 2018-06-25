using ESRI.ArcGIS.Geodatabase;
using ESRI.ArcGIS.Geometry;
using System;
using System.Collections;
using System.IO;


namespace VECTOR_MODEL
{
    class VectorModel
    {
        #region Prepare
        //---------------------------------------------------------------------------------------------
        //path
        string workspacepath = @"C:\Users\B410\Desktop\beforeSchool\optimization\Result3-1-1.5.gdb";
        string eachtinpath = @"C:\Users\B410\Desktop\beforeSchool\optimization\FdISOsurface3-1-1.5\tin_node_";
        string MFilepath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\MFile.txt";
        string MEFilepath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\MEFile.txt";
        string GFilepath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\GFile.txt";
        string GEFilepath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\GEFile.txt";
        string eachgocadpath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\PointsFile";
        string NVXYZpath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\NFile.txt";
        string E1XYZpath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\E1File.txt";
        string E2XYZpath = @"C:\Users\B410\Desktop\beforeSchool\optimization\TXT3-1-1.5\E2File.txt";

        public static int TinCount = 51;                                                                                           //tin网个数，最后一个为初始曲面
        double RasterSize = 5;                                                                                                        //栅格单元边长    
        public static double XScale = 0.002;                                                                                     //
        public static double YScale = 0.002;                                                                                     //箭头放缩
        public static double ZScale = 0.002;                                                                                     //
        public static double NLength_D = 5;                                                                                  //等值面和初始面法矢量长度                                                                                         
        public double NLength_F = 7.5;                                                                                         //集中法矢量长度
        public static double ELengh_D = 100;                                                                                //等值面和初始面曲率矢量长度
        public static double ELengh_F = 150;                                                                                   //集中曲率矢量长度  
        public double ConnectLength = 0.7;                                                                                   //首尾相连步长

        //---------------------------------------------------------------------------------------------
        long TriangleCount;
        long NodeCount;
        IVector3D[] NormalVectorOfNode;
        IVector3D[] NormalVectorOfSurface;
        IPoint[] Centroid;
        ITinAdvanced2 tin;
        //一阶邻域
        public static double NodeRange = 1.2;
        public static double CentriodRange = 1;
        IWorkspace Workspace;
        ITinNode[] DisOrderNode;
        IPoint[] DisOrderPoint;
        double[] KG;
        double[] M;
        double[] k1;
        double[] k2;
        IVector3D[] e1;
        IVector3D[] e2;
        public int D;
        int RowCount;
        double XMin;
        double XMax;
        double YMin;
        double YMax;
        IVector3D[,] NormalVectors;
        //单位主方向
        IVector3D[,] e1s;
        IVector3D[,] e2s;
        IPolyline[,] NormalVector_F;
        IPolyline[,] NormalVector_D;
        IMultiPatch[,] NormalArrow_F;
        IMultiPatch[,] NormalArrow_D;
        IPolyline[,] e1_F;
        IPolyline[,] e2_F;
        IPolyline[,] e1_D;
        IPolyline[,] e2_D;
        IMultiPatch[,] e1Arrow_F;
        IMultiPatch[,] e2Arrow_F;
        IMultiPatch[,] e1Arrow_D;
        IMultiPatch[,] e2Arrow_D;
        double[,] k1_G;
        double[,] k2_G;
        double[,] G_G;
        double[,] M_G;
        double[,] ME_G;
        double[,] GE_G;
        IPoint[,] Points_All;
        //曲率正负判定
        int IsRegular;
        //要素类声明
        //Disperse
        IFields Nodefields_D;
        IFeatureClass NodeFeatureClass_D;
        IFields NormalVectorLinefields_D;
        IFields NormalVectorArrowfields_D;
        IFeatureClass NormalVectorLineFeatureClass_D;
        IFeatureClass NormalVectorArrowFeatureClass_D;
        IFields e1Linefields_D;
        IFields e1Arrowfields_D;
        IFeatureClass e1LineFeatureClass_D;
        IFeatureClass e1ArrowFeatureClass_D;
        IFields e2Linefields_D;
        IFields e2Arrowfields_D;
        IFeatureClass e2LineFeatureClass_D;
        IFeatureClass e2ArrowFeatureClass_D;
        //Foucus
        IFields NormalVectorLinefields_F;
        IFields NormalVectorArrowfields_F;
        IFeatureClass NormalVectorLineFeatureClass_F;
        IFeatureClass NormalVectorArrowFeatureClass_F;
        IFields e1Linefields_F;
        IFields e1Arrowfields_F;
        IFeatureClass e1LineFeatureClass_F;
        IFeatureClass e1ArrowFeatureClass_F;
        IFields e2Linefields_F;
        IFields e2Arrowfields_F;
        IFeatureClass e2LineFeatureClass_F;
        IFeatureClass e2ArrowFeatureClass_F;
        IFields ConnectSLfields;
        IFields ConnectSLArrowfields;
        IFeatureClass ConnectSLFeatureClass;
        IFeatureClass ConnectSLArrowFeatureClass;

        private static object _missing = Type.Missing;

        // 获得TIN数据集
        public void GetTinData(int TinID)
        {
            Type factoryType = Type.GetTypeFromProgID("esriDataSourcesGDB.FileGDBWorkspaceFactory");
            IWorkspaceFactory workspaceFactory = (IWorkspaceFactory)Activator.CreateInstance(factoryType);
            Workspace = workspaceFactory.OpenFromFile(workspacepath, 0);
            // Open an existing TIN, edit the following path as appropriate.
            tin = new TinClass();
            tin.Init(eachtinpath + TinID);
            //得到和tin网直接相关的数据
            TriangleCount = tin.DataTriangleCount;
            NodeCount = tin.NodeCount;
            DisOrderNode = new ITinNode[NodeCount];
            DisOrderPoint = new IPoint[NodeCount];
            //声明预备属性及数组
            if (D == 0)
            {
                IEnvelope Envelope = tin.Extent;
                XMin = Envelope.XMin;
                XMax = Envelope.XMax;
                YMin = Envelope.YMin;
                YMax = Envelope.YMax;
                //列数
                RowCount = (int)((XMax - XMin + 0.5 * RasterSize) / RasterSize) + 1;
                Points_All = new IPoint[NodeCount, TinCount];
                NormalArrow_F = new IMultiPatch[NodeCount, TinCount];
                NormalArrow_D = new IMultiPatch[NodeCount, TinCount];
                NormalVectors = new IVector3D[NodeCount, TinCount];
                //单位主方向
                e1s = new IVector3D[NodeCount, TinCount];
                e2s = new IVector3D[NodeCount, TinCount];
                NormalVector_F = new IPolyline[NodeCount, TinCount];
                NormalVector_D = new IPolyline[NodeCount, TinCount];
                e1_F = new IPolyline[NodeCount, TinCount];
                e2_F = new IPolyline[NodeCount, TinCount];
                e1_D = new IPolyline[NodeCount, TinCount];
                e2_D = new IPolyline[NodeCount, TinCount];
                e1Arrow_F = new IMultiPatch[NodeCount, TinCount];
                e2Arrow_F = new IMultiPatch[NodeCount, TinCount];
                e1Arrow_D = new IMultiPatch[NodeCount, TinCount];
                e2Arrow_D = new IMultiPatch[NodeCount, TinCount];
                k1_G = new double[NodeCount, TinCount];
                k2_G = new double[NodeCount, TinCount];
                G_G = new double[NodeCount, TinCount];
                M_G = new double[NodeCount, TinCount];
                ME_G = new double[NodeCount, TinCount];
                GE_G = new double[NodeCount, TinCount];
            }

            for (int i = 0; i < NodeCount; i++)
            {
                //所有点集
                DisOrderNode[i] = tin.GetNode(i + 1);
                DisOrderPoint[i] = new PointClass();
                DisOrderNode[i].QueryAsPoint(DisOrderPoint[i]);
                Material.MakeZAware(DisOrderPoint[i] as IGeometry);
            }
        }

        ////得到平均tin的曲面要素    
        IPoint[] StartPoints;
        public void GetMeanInfo()
        {
            StartPoints = new IPoint[NodeCount];
            for (int i = 0; i < NodeCount; i++)
            {
                StartPoints[i] = GetOriginal.OPoint[i];
            }
        }
        #endregion Prepare

        #region NormalVector
        //获取三角形的法矢量
        public void GetSurfaceNormalVector()
        {
            ITinAdvanced TinAdvance = tin as ITinAdvanced;
            Centroid = new IPoint[TriangleCount];
            NormalVectorOfSurface = new IVector3D[TriangleCount];
            for (int i = 0; i < TriangleCount; i++)
            {
                Centroid[i] = new PointClass();
                Material.MakeZAware(Centroid[i] as IGeometry);
                NormalVectorOfSurface[i] = new Vector3DClass();
                ITinTriangle TinTriangle = TinAdvance.GetTriangle(i + 1);
                if (TinTriangle.IsEmpty)
                    continue;
                //得到三角形形心
                TinTriangle.QueryCentroid(Centroid[i]);
                //得到基于面的单位法矢量
                TinTriangle.QueryNormal(NormalVectorOfSurface[i]);
                if (NormalVectorOfSurface[i].IsEmpty)
                    continue;
                NormalVectorOfSurface[i].Normalize();
            }
        }

        //获取网格节点p的法矢量
        public void GetNodeNormalVector()
        {
            ITinAdvanced TinAdvance = tin as ITinAdvanced;
            NormalVectorOfNode = new IVector3D[NodeCount];
            //对每个节点进行循环
            for (int i = 4; i < NodeCount; i++)
            {
                NormalVectorOfNode[i] = new Vector3DClass();
                ITinTriangleArray IncidentTriangles = Material.GetIncidentTriangles(DisOrderNode[i], CentriodRange);
                int IncidentTriangleCount = IncidentTriangles.Count;
                //对p节点的一阶邻域节点循环
                for (int j = 0; j < IncidentTriangleCount; j++)
                {
                    IPoint IncidentCentroid = new PointClass();
                    Material.MakeZAware(IncidentCentroid as IGeometry);
                    IVector3D IncidentUnitNormalVectorOfSurface = new Vector3DClass();
                    ITinTriangle IncidentTinTriangle = IncidentTriangles.get_Element(j);
                    if (IncidentTinTriangle.IsEmpty)
                        continue;
                    IncidentTinTriangle.QueryNormal(IncidentUnitNormalVectorOfSurface);
                    if (IncidentUnitNormalVectorOfSurface.IsEmpty)
                        continue;
                    //得到p节点相邻三角形的形心和单位法矢量,得到节点到形心的距离
                    IncidentTinTriangle.QueryCentroid(IncidentCentroid);
                    double NCdistance = System.Math.Sqrt((IncidentCentroid.X - DisOrderNode[i].X) * (IncidentCentroid.X - DisOrderNode[i].X) +
                    (IncidentCentroid.Y - DisOrderNode[i].Y) * (IncidentCentroid.Y - DisOrderNode[i].Y) + (IncidentCentroid.Z - DisOrderNode[i].Z) * (IncidentCentroid.Z - DisOrderNode[i].Z));
                    //得到节点的法矢量分量,合成法矢量
                    IncidentUnitNormalVectorOfSurface.Magnitude = (1 / NCdistance);
                    if ((NormalVectorOfNode[i] == null) || (NormalVectorOfNode[i].IsEmpty))
                    {
                        NormalVectorOfNode[i] = IncidentUnitNormalVectorOfSurface;
                    }
                    else
                    {
                        NormalVectorOfNode[i].ConstructAddVector(NormalVectorOfNode[i], IncidentUnitNormalVectorOfSurface);
                    }
                }
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                NormalVectorOfNode[i].Normalize();
            }
        }
        #endregion  NormalVector

        #region DifferentialGeometry
        public void GetNodeCurvature()
        {
            //将向量 pi−p 投影到 p 点所在的切平面,得到切平面的切向单位向量ti
            ITinAdvanced TinAdvance = tin as ITinAdvanced;
            IVector3D[] tid = new IVector3D[NodeCount];
            double[] a = new double[NodeCount];
            double[] b = new double[NodeCount];
            double[] c = new double[NodeCount];
            IVector3D[] E1 = new IVector3D[NodeCount];
            IVector3D[] E2 = new IVector3D[NodeCount];
            double[] a11 = new double[NodeCount];
            double[] a12 = new double[NodeCount];
            double[] a21 = new double[NodeCount];
            double[] a22 = new double[NodeCount];
            double[] a13 = new double[NodeCount];
            double[] a23 = new double[NodeCount];
            double[] O = new double[NodeCount];
            KG = new double[NodeCount];
            M = new double[NodeCount];
            k1 = new double[NodeCount];
            k2 = new double[NodeCount];
            e1 = new IVector3D[NodeCount];
            e2 = new IVector3D[NodeCount];

            //对每个节点进行循环
            for (int i = 4; i < NodeCount; i++)
            {
                ITinNodeArray AdjacentNodesArray = Material.GetIncidentNodes(DisOrderNode[i], NodeRange);
                int AdjacentNodesCount = AdjacentNodesArray.Count;

                //沿着edge的切向向量ti
                IVector3D[] ti = new IVector3D[AdjacentNodesCount];
                //沿ti方向的法曲率Knti
                double[] Knti = new double[AdjacentNodesCount];

                //对p节点的一阶邻域循环得到edge曲率
                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    //求ti
                    ITinNode AdjacentNode = AdjacentNodesArray.get_Element(j);
                    //转换Node为Point
                    IPoint AdjacentPoint = new PointClass();
                    AdjacentNode.QueryAsPoint(AdjacentPoint);
                    Material.MakeZAware(AdjacentPoint as IGeometry);
                    //(pi-p)
                    IVector3D EdgeVectorsOfNode = Material.CreateVector3DTwoPoints(AdjacentPoint, DisOrderPoint[i]);
                    //<pi-p,N>N
                    IVector3D ProjectNormalVector = new Vector3DClass();
                    Material.Clone(NormalVectorOfNode[i], ProjectNormalVector);
                    IVector3D NormalVectorOfNode2 = new Vector3DClass();
                    Material.Clone(NormalVectorOfNode[i], NormalVectorOfNode2);
                    if (ProjectNormalVector.IsEmpty)
                        continue;
                    if (NormalVectorOfNode[i].IsEmpty)
                        continue;
                    ProjectNormalVector.Magnitude = EdgeVectorsOfNode.DotProduct(NormalVectorOfNode2);
                    //ti=(pi-p)-<pi-p,N>N
                    ti[j] = new Vector3DClass();
                    ti[j].ConstructSubtractVector(EdgeVectorsOfNode, ProjectNormalVector);
                    ti[j].Normalize();
                    //Ni-N
                    IVector3D NiN = new Vector3DClass();
                    if (NormalVectorOfNode[(AdjacentNode.Index) - 1] == null)
                        continue;
                    if (NormalVectorOfNode[(AdjacentNode.Index) - 1].IsEmpty)
                        continue;
                    if (NormalVectorOfNode[i].IsEmpty)
                        continue;
                    NiN.ConstructSubtractVector(NormalVectorOfNode[(AdjacentNode.Index) - 1], NormalVectorOfNode[i]);
                    //Kn(ti)
                    Knti[j] = -(EdgeVectorsOfNode.DotProduct(NiN)) / (EdgeVectorsOfNode.DotProduct(EdgeVectorsOfNode));
                }
                //得到最大的Knti曲率Kntid
                Random r = new Random();
                int e = r.Next(AdjacentNodesCount);
                double Kntid = Knti[e];
                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    if (Knti[j] > Kntid)
                    {
                        Kntid = Knti[j];
                    }
                    a[i] = Kntid;
                }
                //得到节点p的最大Kntid所对应的ti，即tid
                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    if (Knti[j] == Kntid)
                    {
                        tid[i] = ti[j];
                        break;
                    }
                }
                //得到坐标轴e1,e2,夹角cos(&i),sin(&i),以及系数aij
                E1[i] = tid[i];
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                E2[i] = new Vector3DClass();
                E2[i].ConstructCrossProduct(E1[i], NormalVectorOfNode[i]);
                E2[i].Normalize();
                //计算各项系数
                double[] A11 = new double[AdjacentNodesCount];
                double[] A12 = new double[AdjacentNodesCount];
                double[] A21 = new double[AdjacentNodesCount];
                double[] A22 = new double[AdjacentNodesCount];
                double[] A13 = new double[AdjacentNodesCount];
                double[] A23 = new double[AdjacentNodesCount];

                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    double cos = ti[j].DotProduct(E1[i]);
                    double sin = ti[j].DotProduct(E2[i]);
                    A11[j] = cos * cos * sin * sin;
                    A12[j] = cos * sin * sin * sin;
                    A21[j] = A12[j];
                    A22[j] = sin * sin * sin * sin;
                    A13[j] = (Knti[j] - a[i] * cos * cos) * cos * sin;
                    A23[j] = (Knti[j] - a[i] * cos * cos) * sin * sin;
                }

                a11[i] = a12[i] = a21[i] = a22[i] = a13[i] = a23[i] = 0;
                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    a11[i] += A11[j];
                    a12[i] += A12[j];
                    a21[i] += A21[j];
                    a22[i] += A22[j];
                    a13[i] += A13[j];
                    a23[i] += A23[j];
                }

                //得到系数b，c
                b[i] = (a13[i] * a22[i] - a23[i] * a12[i]) / (a11[i] * a22[i] - a12[i] * a12[i]);
                c[i] = (a11[i] * a23[i] - a12[i] * a13[i]) / (a11[i] * a22[i] - a12[i] * a12[i]);

                //高斯曲率，平均曲率，主曲率
                KG[i] = a[i] * c[i] - (b[i] * b[i] / 4);
                M[i] = (a[i] + c[i]) / 2;
                k1[i] = M[i] + Math.Sqrt(M[i] * M[i] - KG[i]);
                k2[i] = M[i] - Math.Sqrt(M[i] * M[i] - KG[i]);

                //主方向E1，E2，其中O为计算角度
                O[i] = 0.5 * Math.Asin(b[i] / (k2[i] - k1[i]));
                IVector3D tempE1 = new Vector3DClass();
                IVector3D tranE1 = new Vector3DClass();
                IVector3D momE1 = new Vector3DClass();
                IVector3D briE1 = new Vector3DClass();
                IVector3D tempE2 = new Vector3DClass();
                IVector3D tranE2 = new Vector3DClass();
                IVector3D momE2 = new Vector3DClass();
                IVector3D briE2 = new Vector3DClass();
                if (E1[i].IsEmpty)
                    continue;
                if (E2[i].IsEmpty)
                    continue;
                Material.Clone(E1[i], tempE1);
                Material.Clone(E2[i], tempE2);
                Material.Clone(E1[i], tranE1);
                Material.Clone(E2[i], tranE2);
                Material.Clone(E1[i], momE1);
                Material.Clone(E2[i], momE2);
                Material.Clone(E1[i], briE1);
                Material.Clone(E2[i], briE2);
                //e1
                tempE1.Magnitude = (momE1.Magnitude) * (Math.Cos(O[i]));
                tempE2.Magnitude = (momE2.Magnitude) * (Math.Sin(O[i]));
                e1[i] = new Vector3DClass();
                if (tempE1.IsEmpty)
                    continue;
                if (tempE2.IsEmpty)
                    continue;
                e1[i].ConstructAddVector(tempE1, tempE2);
                //e2
                tranE1.Magnitude = ((briE1.Magnitude) * (Math.Sin(O[i])));
                tranE2.Magnitude = ((briE2.Magnitude) * (Math.Cos(O[i])));
                e2[i] = new Vector3DClass();
                if (tranE1.IsEmpty)
                    continue;
                if (tranE1.IsEmpty)
                    continue;
                e2[i].ConstructSubtractVector(tranE2, tranE1);
            }
        }
        #endregion DifferentialGeometry

        #region CreateFeature
        //先建全局矩阵然后在最后保存
        #region Prepare
        //要素类的名称      
        string NodeClass_D = "Node_D";
        string NormalVectorLineClass_D = "NormalVector_D";
        string NormalVectorArrowClass_D = "NormalVectorArrow_D";
        string e1LineClass_D = "e1_D";
        string e1ArrowClass_D = "e1Arrow_D";
        string e2LineClass_D = "e2_D";
        string e2ArrowClass_D = "e2Arrow_D";
        string NormalVectorLineClass_F = "NormalVector_F";
        string NormalVectorArrowClass_F = "NormalVectorArrow1_F";
        string e1LineClass_F = "e1_F";
        string e1ArrowClass_F = "e1Arrow_F";
        string e2LineClass_F = "e2_F";
        string e2ArrowClass_F = "e2Arrow_F";
        string ConnectSLClass = "ConnectSL";
        string ConnectSLArrowClass = "ConnectSLArrow";

        //要素类
        public IFeatureClass CreateStandaloneFeatureClass(IWorkspace workspace, String featureClassName,
        IFields fieldsCollection)
        {
            IFeatureWorkspace featureWorkspace = (IFeatureWorkspace)workspace;
            IFeatureClassDescription fcDesc = new FeatureClassDescriptionClass();
            IObjectClassDescription ocDesc = (IObjectClassDescription)fcDesc;

            // Use IFieldChecker to create a validated fields collection.
            IFieldChecker fieldChecker = new FieldCheckerClass();
            IEnumFieldError enumFieldError = null;
            IFields validatedFields = null;
            fieldChecker.ValidateWorkspace = workspace;
            fieldChecker.Validate(fieldsCollection, out enumFieldError, out validatedFields);

            // The enumFieldError enumerator can be inspected at this point to determine 
            // which fields were modified during validation.
            IFeatureClass featureClass = featureWorkspace.CreateFeatureClass(featureClassName, validatedFields,
            ocDesc.InstanceCLSID, ocDesc.ClassExtensionCLSID, esriFeatureType.esriFTSimple, fcDesc.ShapeFieldName, "");
            return featureClass;
        }

        //字段集
        private IFields CreateFieldsCollectionForFeatueClass(ISpatialReference spatialReference,
            esriGeometryType geometryType, string className)
        {
            //Use the feature class description to return the required fields in a fields collection.
            IFeatureClassDescription fcDesc = new FeatureClassDescriptionClass();
            IObjectClassDescription ocDesc = (IObjectClassDescription)fcDesc;

            //Create the fields using the required fields method.
            IFields fields = ocDesc.RequiredFields;

            //Locate the shape field with the name from the feature class description.
            int shapeFieldIndex = fields.FindField(fcDesc.ShapeFieldName);
            IField shapeField = fields.get_Field(shapeFieldIndex);

            //Modify the GeometryDef object before using the fields collection to create a
            //feature class.
            IGeometryDef geometryDef = shapeField.GeometryDef;
            IGeometryDefEdit geometryDefEdit = (IGeometryDefEdit)geometryDef;

            //Alter the feature class geometry type to the type we need.
            geometryDefEdit.GeometryType_2 = geometryType;
            geometryDefEdit.HasZ_2 = true;
            geometryDefEdit.HasM_2 = true;
            geometryDefEdit.GridCount_2 = 1;

            //Set the first grid size to zero and allow ArcGIS to determine a valid grid size.
            geometryDefEdit.set_GridSize(0, 0);
            geometryDefEdit.SpatialReference_2 = spatialReference;

            //If it is a point collection, add Field
            if (geometryType == esriGeometryType.esriGeometryPoint)
            {
                //Create a user-defined double field for "Gauss"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "Gauss";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "GaussCurvature";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
                //Create a user-defined double field for "H"
                IFieldsEdit fieldsEdit2 = (IFieldsEdit)fields;
                IField incomeField2 = new FieldClass();
                IFieldEdit incomeFieldEdit2 = (IFieldEdit)incomeField2;
                incomeFieldEdit2.AliasName_2 = "H";
                incomeFieldEdit2.Editable_2 = true;
                incomeFieldEdit2.IsNullable_2 = false;
                incomeFieldEdit2.Name_2 = "MeanCurvature";
                incomeFieldEdit2.Precision_2 = 2;
                incomeFieldEdit2.Scale_2 = 5;
                incomeFieldEdit2.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit2.AddField(incomeField2);
                //Create a user-defined double field for "k1"
                IFieldsEdit fieldsEdit3 = (IFieldsEdit)fields;
                IField incomeField3 = new FieldClass();
                IFieldEdit incomeFieldEdit3 = (IFieldEdit)incomeField3;
                incomeFieldEdit3.AliasName_2 = "k1";
                incomeFieldEdit3.Editable_2 = true;
                incomeFieldEdit3.IsNullable_2 = false;
                incomeFieldEdit3.Name_2 = "BigCurvature";
                incomeFieldEdit3.Precision_2 = 2;
                incomeFieldEdit3.Scale_2 = 5;
                incomeFieldEdit3.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit3.AddField(incomeField3);
                //Create a user-defined double field for "k2"
                IFieldsEdit fieldsEdit4 = (IFieldsEdit)fields;
                IField incomeField4 = new FieldClass();
                IFieldEdit incomeFieldEdit4 = (IFieldEdit)incomeField4;
                incomeFieldEdit4.AliasName_2 = "k2";
                incomeFieldEdit4.Editable_2 = true;
                incomeFieldEdit4.IsNullable_2 = false;
                incomeFieldEdit4.Name_2 = "SmallCurvature";
                incomeFieldEdit4.Precision_2 = 2;
                incomeFieldEdit4.Scale_2 = 5;
                incomeFieldEdit4.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit4.AddField(incomeField4);
                //Create a user-defined double field for "X"
                IFieldsEdit fieldsEdit5 = (IFieldsEdit)fields;
                IField incomeField5 = new FieldClass();
                IFieldEdit incomeFieldEdit5 = (IFieldEdit)incomeField5;
                incomeFieldEdit5.AliasName_2 = "X";
                incomeFieldEdit5.Editable_2 = true;
                incomeFieldEdit5.IsNullable_2 = false;
                incomeFieldEdit5.Name_2 = "X";
                incomeFieldEdit5.Precision_2 = 2;
                incomeFieldEdit5.Scale_2 = 5;
                incomeFieldEdit5.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit5.AddField(incomeField5);
                //Create a user-defined double field for "Y"
                IFieldsEdit fieldsEdit6 = (IFieldsEdit)fields;
                IField incomeField6 = new FieldClass();
                IFieldEdit incomeFieldEdit6 = (IFieldEdit)incomeField6;
                incomeFieldEdit6.AliasName_2 = "Y";
                incomeFieldEdit6.Editable_2 = true;
                incomeFieldEdit6.IsNullable_2 = false;
                incomeFieldEdit6.Name_2 = "Y";
                incomeFieldEdit6.Precision_2 = 2;
                incomeFieldEdit6.Scale_2 = 5;
                incomeFieldEdit6.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit6.AddField(incomeField6);
                //Create a user-defined double field for "Z"
                IFieldsEdit fieldsEdit7 = (IFieldsEdit)fields;
                IField incomeField7 = new FieldClass();
                IFieldEdit incomeFieldEdit7 = (IFieldEdit)incomeField7;
                incomeFieldEdit7.AliasName_2 = "Z";
                incomeFieldEdit7.Editable_2 = true;
                incomeFieldEdit7.IsNullable_2 = false;
                incomeFieldEdit7.Name_2 = "Z";
                incomeFieldEdit7.Precision_2 = 2;
                incomeFieldEdit7.Scale_2 = 5;
                incomeFieldEdit7.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit7.AddField(incomeField7);
                //Create a user-defined double field for "Probability"
                IFieldsEdit fieldsEdit148 = (IFieldsEdit)fields;
                IField incomeField148 = new FieldClass();
                IFieldEdit incomeFieldEdit148 = (IFieldEdit)incomeField148;
                incomeFieldEdit148.AliasName_2 = "Probability";
                incomeFieldEdit148.Editable_2 = true;
                incomeFieldEdit148.IsNullable_2 = false;
                incomeFieldEdit148.Name_2 = "Probability";
                incomeFieldEdit148.Precision_2 = 2;
                incomeFieldEdit148.Scale_2 = 5;
                incomeFieldEdit148.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit148.AddField(incomeField148);
                //Create a user-defined double field for "HEntropy"
                IFieldsEdit fieldsEdit019 = (IFieldsEdit)fields;
                IField incomeField019 = new FieldClass();
                IFieldEdit incomeFieldEdit019 = (IFieldEdit)incomeField019;
                incomeFieldEdit019.AliasName_2 = "MEntropy";
                incomeFieldEdit019.Editable_2 = true;
                incomeFieldEdit019.IsNullable_2 = false;
                incomeFieldEdit019.Name_2 = "MEntropy";
                incomeFieldEdit019.Precision_2 = 2;
                incomeFieldEdit019.Scale_2 = 5;
                incomeFieldEdit019.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit019.AddField(incomeField019);
                //Create a user-defined double field for "GEntropy"
                IFieldsEdit fieldsEdit0119 = (IFieldsEdit)fields;
                IField incomeField0119 = new FieldClass();
                IFieldEdit incomeFieldEdit0119 = (IFieldEdit)incomeField0119;
                incomeFieldEdit0119.AliasName_2 = "GEntropy";
                incomeFieldEdit0119.Editable_2 = true;
                incomeFieldEdit0119.IsNullable_2 = false;
                incomeFieldEdit0119.Name_2 = "GEntropy";
                incomeFieldEdit0119.Precision_2 = 2;
                incomeFieldEdit0119.Scale_2 = 5;
                incomeFieldEdit0119.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit0119.AddField(incomeField0119);
            }

            //If it is a NormalVector, add Field
            if ((className == NormalVectorLineClass_F) || (className == NormalVectorLineClass_D))
            {
                //Create a user-defined double field for "Magnitude"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "Magnitude";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "Magnitude";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
                //Create a user-defined double field for "Probability"
                IFieldsEdit fieldsEdit148 = (IFieldsEdit)fields;
                IField incomeField148 = new FieldClass();
                IFieldEdit incomeFieldEdit148 = (IFieldEdit)incomeField148;
                incomeFieldEdit148.AliasName_2 = "Probability";
                incomeFieldEdit148.Editable_2 = true;
                incomeFieldEdit148.IsNullable_2 = false;
                incomeFieldEdit148.Name_2 = "Probability";
                incomeFieldEdit148.Precision_2 = 2;
                incomeFieldEdit148.Scale_2 = 5;
                incomeFieldEdit148.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit148.AddField(incomeField148);
            }
            //If it is a ConnectSL, add Field
            if (className == ConnectSLClass)
            {
                //Create a user-defined double field for "Probability"
                IFieldsEdit fieldsEdit148 = (IFieldsEdit)fields;
                IField incomeField148 = new FieldClass();
                IFieldEdit incomeFieldEdit148 = (IFieldEdit)incomeField148;
                incomeFieldEdit148.AliasName_2 = "Probability";
                incomeFieldEdit148.Editable_2 = true;
                incomeFieldEdit148.IsNullable_2 = false;
                incomeFieldEdit148.Name_2 = "Probability";
                incomeFieldEdit148.Precision_2 = 2;
                incomeFieldEdit148.Scale_2 = 5;
                incomeFieldEdit148.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit148.AddField(incomeField148);
            }

            if ((className == e1LineClass_F) || (className == e1LineClass_D))
            {
                //Create a user-defined double field for "k1"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "k1";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "k1";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
            }
            if ((className == e2LineClass_F) || (className == e2LineClass_D))
            {
                //Create a user-defined double field for "k2"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "k2";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "k2";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
            }

            if ((className == e1LineClass_F) || (className == e2LineClass_F) || (className == e1LineClass_D) || (className == e2LineClass_D))
            {
                //Create a user-defined double field for "Probability"
                IFieldsEdit fieldsEdit148 = (IFieldsEdit)fields;
                IField incomeField148 = new FieldClass();
                IFieldEdit incomeFieldEdit148 = (IFieldEdit)incomeField148;
                incomeFieldEdit148.AliasName_2 = "Probability";
                incomeFieldEdit148.Editable_2 = true;
                incomeFieldEdit148.IsNullable_2 = false;
                incomeFieldEdit148.Name_2 = "Probability";
                incomeFieldEdit148.Precision_2 = 2;
                incomeFieldEdit148.Scale_2 = 5;
                incomeFieldEdit148.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit148.AddField(incomeField148);
            }

            if ((className == e1ArrowClass_F) || (className == e2ArrowClass_F) || (className == e1ArrowClass_D) || (className == e2ArrowClass_D))
            {
                //Create a user-defined double field for "IsRegular"
                IFieldsEdit fieldsEdit3 = (IFieldsEdit)fields;
                IField incomeField3 = new FieldClass();
                IFieldEdit incomeFieldEdit3 = (IFieldEdit)incomeField3;
                incomeFieldEdit3.AliasName_2 = "IsRegular";
                incomeFieldEdit3.Editable_2 = true;
                incomeFieldEdit3.IsNullable_2 = false;
                incomeFieldEdit3.Name_2 = "IsRegular";
                incomeFieldEdit3.Precision_2 = 2;
                incomeFieldEdit3.Scale_2 = 5;
                incomeFieldEdit3.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit3.AddField(incomeField3);
            }
            return fields;
        }

        IGeoDataset geoDataset;
        ISpatialReference spatialReference;
        public void CreateFeatureClass()
        {
            //传递坐标系
            geoDataset = tin as IGeoDataset;
            spatialReference = geoDataset.SpatialReference;
            //创建要素类
            //Disperse
            Nodefields_D = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPoint, NodeClass_D);
            NodeFeatureClass_D = CreateStandaloneFeatureClass(Workspace, NodeClass_D, Nodefields_D);

            NormalVectorLinefields_D = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline, NormalVectorLineClass_D);
            NormalVectorArrowfields_D = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch, NormalVectorArrowClass_D);
            NormalVectorLineFeatureClass_D = CreateStandaloneFeatureClass(Workspace, NormalVectorLineClass_D, NormalVectorLinefields_D);
            NormalVectorArrowFeatureClass_D = CreateStandaloneFeatureClass(Workspace, NormalVectorArrowClass_D, NormalVectorArrowfields_D);

            e1Linefields_D = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline, e1LineClass_D);
            e1Arrowfields_D = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch, e1ArrowClass_D);
            e1LineFeatureClass_D = CreateStandaloneFeatureClass(Workspace, e1LineClass_D, e1Linefields_D);
            e1ArrowFeatureClass_D = CreateStandaloneFeatureClass(Workspace, e1ArrowClass_D, e1Arrowfields_D);

            e2Linefields_D = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline, e2LineClass_D);
            e2Arrowfields_D = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch, e2ArrowClass_D);
            e2LineFeatureClass_D = CreateStandaloneFeatureClass(Workspace, e2LineClass_D, e2Linefields_D);
            e2ArrowFeatureClass_D = CreateStandaloneFeatureClass(Workspace, e2ArrowClass_D, e2Arrowfields_D);

            //Focus
            NormalVectorLinefields_F = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline, NormalVectorLineClass_F);
            NormalVectorArrowfields_F = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch, NormalVectorArrowClass_F);
            NormalVectorLineFeatureClass_F = CreateStandaloneFeatureClass(Workspace, NormalVectorLineClass_F, NormalVectorLinefields_F);
            NormalVectorArrowFeatureClass_F = CreateStandaloneFeatureClass(Workspace, NormalVectorArrowClass_F, NormalVectorArrowfields_F);

            e1Linefields_F = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline, e1LineClass_F);
            e1Arrowfields_F = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch, e1ArrowClass_F);
            e1LineFeatureClass_F = CreateStandaloneFeatureClass(Workspace, e1LineClass_F, e1Linefields_F);
            e1ArrowFeatureClass_F = CreateStandaloneFeatureClass(Workspace, e1ArrowClass_F, e1Arrowfields_F);

            e2Linefields_F = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline, e2LineClass_F);
            e2Arrowfields_F = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch, e2ArrowClass_F);
            e2LineFeatureClass_F = CreateStandaloneFeatureClass(Workspace, e2LineClass_F, e2Linefields_F);
            e2ArrowFeatureClass_F = CreateStandaloneFeatureClass(Workspace, e2ArrowClass_F, e2Arrowfields_F);

            //SL          
            ConnectSLfields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline, ConnectSLClass);
            ConnectSLArrowfields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch, ConnectSLArrowClass);
            ConnectSLFeatureClass = CreateStandaloneFeatureClass(Workspace, ConnectSLClass, ConnectSLfields);
            ConnectSLArrowFeatureClass = CreateStandaloneFeatureClass(Workspace, ConnectSLArrowClass, ConnectSLArrowfields);
        }

        public void CreateFeature()
        {
            //移动放缩矢量箭头
            IVector3D VerticalVector = Material.ConstructVector3D(0, 0, 1);
            IPoint originPoint = Material.ConstructPoint3D(0, 0, 0);
            Material.MakeZAware(originPoint as IGeometry);
        #endregion Prepare

            #region Node
            //构建顶点及其层序数组
            for (int i = 4; i < NodeCount; i++)
            {
                //存入数组待比较
                int v = NodeLocation(DisOrderNode[i]); ;
                Points_All[v, D] = DisOrderPoint[i];
                Material.MakeZAware(Points_All[v, D]);
                //G和H
                G_G[v, D] = KG[i];
                M_G[v, D] = M[i];
            }

            //熵
            //为H，G，Node排序
            double[] orderH = new double[NodeCount];
            double[] orderG = new double[NodeCount];
            ITinNode[] orderNode = new ITinNode[NodeCount];
            for (int i = 4; i < NodeCount; i++)
            {
                int v = NodeLocation(DisOrderNode[i]);
                orderH[v] = M[i];
                orderG[v] = KG[i];
                orderNode[v] = DisOrderNode[i];
            }

            for (int i = 0; i < NodeCount - 4; i++)
            {
                ME_G[i, D] = Material.Entropy(orderNode[i], orderH, tin);
                GE_G[i, D] = Material.Entropy(orderNode[i], orderG, tin);
            }
            #endregion Node

            #region NodeNormalVector
            //创建基于顶点的法矢量Focus
            for (int i = 4; i < NodeCount; i++)
            {
                if (DisOrderPoint[i].IsEmpty)
                    continue;
                //更新其ID
                int v = NodeLocation(DisOrderNode[i]);
                IRay ray = new RayClass();
                ray.Origin = StartPoints[v];
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                ray.Vector = NormalVectorOfNode[i];
                //设置矢量可视化长度
                IPoint endPoint = ray.GetPointAtDistance(NLength_F);                                              //cd                     
                Material.MakeZAware(endPoint as IGeometry);
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                //画出矢量线
                VectorPointCollection.AddPoint(StartPoints[v], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);

                //为向量末端添加箭头符号
                IGeometry Arrow1 = Material.GetArrow();
                double Inclination1 = NormalVectorOfNode[i].Inclination;
                double degreesOfRotation1 = Math.Acos((NormalVectorOfNode[i].DotProduct(VerticalVector)) / ((NormalVectorOfNode[i].Magnitude) * (VerticalVector.Magnitude)));                  //
                ITransform3D transform3D1 = Arrow1 as ITransform3D;
                //放缩
                transform3D1.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation1 != 0)
                {
                    double angleOfRotationInRadians1 = degreesOfRotation1;
                    IVector3D axisOfRotationVector3D1 = new Vector3DClass();
                    axisOfRotationVector3D1.ConstructCrossProduct(VerticalVector, NormalVectorOfNode[i]);
                    transform3D1.RotateVector3D(axisOfRotationVector3D1, angleOfRotationInRadians1);
                }
                //平移
                if (endPoint.IsEmpty)
                    continue;
                transform3D1.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //全局矩阵
                NormalVector_F[v, D] = VectorPointCollection as IPolyline;
                NormalArrow_F[v, D] = Arrow1 as IMultiPatch;
            }

            //创建基于顶点的法矢量Disperse
            for (int i = 4; i < NodeCount; i++)
            {
                if (DisOrderPoint[i].IsEmpty)
                    continue;
                //更新其ID
                int v = NodeLocation(DisOrderNode[i]);
                if (DisOrderPoint[i].IsEmpty)
                    continue;
                IRay ray = new RayClass();
                ray.Origin = DisOrderPoint[i];
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                ray.Vector = NormalVectorOfNode[i];
                //存入法矢量云数组
                NormalVectors[v, D] = NormalVectorOfNode[i];
                //设置矢量可视化长度
                IPoint endPoint = ray.GetPointAtDistance(NLength_D);                                                        //
                IPoint SLendPoint = ray.GetPointAtDistance(NLength_D);
                //           
                Material.MakeZAware(endPoint as IGeometry);
                Material.MakeZAware(SLendPoint);
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                //画出矢量线
                VectorPointCollection.AddPoint(DisOrderPoint[i], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                double Inclination = NormalVectorOfNode[i].Inclination;
                double degreesOfRotation = Math.Acos((NormalVectorOfNode[i].DotProduct(VerticalVector)) / ((NormalVectorOfNode[i].Magnitude) * (VerticalVector.Magnitude)));                  //
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, NormalVectorOfNode[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                if (endPoint.IsEmpty)
                    continue;
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //全局矩阵
                NormalVector_D[v, D] = VectorPointCollection as IPolyline;
                NormalArrow_D[v, D] = Arrow as IMultiPatch;
            }
            #endregion NodeNormalVector

            #region PrincipalDirection
            //创建基于顶点的主方向矢量线
            //e1 Focus
            for (int i = 4; i < NodeCount; i++)
            {
                if (DisOrderPoint[i].IsEmpty)
                    continue;
                int v = NodeLocation(DisOrderNode[i]);
                IRay ray = new RayClass();
                ray.Origin = StartPoints[v];
                if (e1[i] == null)
                    continue;
                ray.Vector = e1[i];
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //设置矢量可视化长度
                endPoint = ray.GetPointAtDistance(Math.Abs(ELengh_F * k1[i]));
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(StartPoints[v], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                if (e1[i].IsEmpty)
                    continue;
                double Inclination = e1[i].Inclination;
                double degreesOfRotation = Math.Acos((e1[i].DotProduct(VerticalVector)) / ((e1[i].Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, e1[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //存入全局矩阵
                if (VectorPointCollection == null)
                    continue;
                e1_F[v, D] = VectorPointCollection as IPolyline;
                k1_G[v, D] = k1[i];
                e1Arrow_F[v, D] = Arrow as IMultiPatch;
                //e1s
                e1s[v, D] = e1[i];
            }

            //e1 Disperse
            for (int i = 4; i < NodeCount; i++)
            {
                if (DisOrderPoint[i].IsEmpty)
                    continue;
                int v = NodeLocation(DisOrderNode[i]);
                IRay ray = new RayClass();
                ray.Origin = DisOrderPoint[i];
                if (e1[i] == null)
                    continue;
                ray.Vector = e1[i];
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //设置矢量可视化长度
                endPoint = ray.GetPointAtDistance(Math.Abs(ELengh_D * k1[i]));                                                                         //
                IPoint SLendPoint = new PointClass();
                SLendPoint = ray.GetPointAtDistance(Math.Abs(ELengh_D * k1[i]));                                                                       //
                Material.MakeZAware(SLendPoint);
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(DisOrderPoint[i], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                if (e1[i].IsEmpty)
                    continue;
                double Inclination = e1[i].Inclination;
                double degreesOfRotation = Math.Acos((e1[i].DotProduct(VerticalVector)) / ((e1[i].Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, e1[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //存入全局矩阵
                if (VectorPointCollection == null)
                    continue;
                e1_D[v, D] = VectorPointCollection as IPolyline;
                e1Arrow_D[v, D] = Arrow as IMultiPatch;
            }

            //e2 Focus
            for (int i = 4; i < NodeCount; i++)
            {
                if (DisOrderPoint[i].IsEmpty)
                    continue;
                int v = NodeLocation(DisOrderNode[i]);
                IRay ray = new RayClass();
                ray.Origin = StartPoints[v];
                if (e2[i] == null)
                    continue;
                ray.Vector = e2[i];
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //设置矢量可视化长度
                endPoint = ray.GetPointAtDistance((Math.Abs(ELengh_F * k2[i])));
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(StartPoints[v], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                if (e2[i] == null)
                    continue;
                double Inclination = e2[i].Inclination;
                double degreesOfRotation = Math.Acos((e2[i].DotProduct(VerticalVector)) / ((e2[i].Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, e2[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                if (endPoint.IsEmpty)
                    continue;
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);

                //存入数组
                if (VectorPointCollection == null)
                    continue;
                e2_F[v, D] = VectorPointCollection as IPolyline;
                k2_G[v, D] = k2[i];
                e2Arrow_F[v, D] = Arrow as IMultiPatch;
                //e1s
                e2s[v, D] = e2[i];
            }

            //e2 Disperse
            for (int i = 4; i < NodeCount; i++)
            {
                if (DisOrderPoint[i].IsEmpty)
                    continue;
                int v = NodeLocation(DisOrderNode[i]);
                IRay ray = new RayClass();
                ray.Origin = DisOrderPoint[i];
                if (e2[i] == null)
                    continue;
                ray.Vector = e2[i];
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //设置矢量可视化长度
                endPoint = ray.GetPointAtDistance((Math.Abs(ELengh_D * k2[i])));                                                                                              //            
                IPoint SLendPoint = new PointClass();
                SLendPoint = ray.GetPointAtDistance(Math.Abs(ELengh_D * k2[i]));                                                                                               //
                Material.MakeZAware(SLendPoint);
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(DisOrderPoint[i], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                if (e2[i] == null)
                    continue;
                double Inclination = e2[i].Inclination;
                double degreesOfRotation = Math.Acos((e2[i].DotProduct(VerticalVector)) / ((e2[i].Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, e2[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                if (endPoint.IsEmpty)
                    continue;
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //存入数组
                if (VectorPointCollection == null)
                    continue;
                e2_D[v, D] = VectorPointCollection as IPolyline;
                e2Arrow_D[v, D] = Arrow as IMultiPatch;
            }
            #endregion PrincipalDirection

            #region Creat Global Feature
            if (D == TinCount - 1)
            {
                //创建集中点并保存其各种属性
                for (int i = 0; i < NodeCount; i++)
                {
                    for (int j = 0; j < TinCount - 1; j++)
                    {
                        if ((Points_All[i, j] == null) || (Points_All[i, j].IsEmpty))
                        {
                            continue;
                        }
                        IFeature feature1 = NodeFeatureClass_D.CreateFeature();
                        feature1.Shape = Points_All[i, j];
                        int index11 = feature1.Fields.FindField("X");
                        feature1.set_Value(index11, Points_All[i, j].X);
                        int index12 = feature1.Fields.FindField("Y");
                        feature1.set_Value(index12, Points_All[i, j].Y);
                        int index13 = feature1.Fields.FindField("Z");
                        feature1.set_Value(index13, Points_All[i, j].Z);
                        int index3 = feature1.Fields.FindField("BigCurvature");
                        feature1.set_Value(index3, k1_G[i, j]);
                        int index4 = feature1.Fields.FindField("SmallCurvature");
                        feature1.set_Value(index4, k2_G[i, j]);
                        int index1 = feature1.Fields.FindField("GaussCurvature");
                        feature1.set_Value(index1, G_G[i, j]);
                        int index2 = feature1.Fields.FindField("MeanCurvature");
                        feature1.set_Value(index2, M_G[i, j]);
                        int index167 = feature1.Fields.FindField("Probability");
                        feature1.set_Value(index167, 0.95 - (0.9 / 49) * j);
                        int index267 = feature1.Fields.FindField("MEntropy");
                        feature1.set_Value(index267, ME_G[i, j]);
                        int index1267 = feature1.Fields.FindField("GEntropy");
                        feature1.set_Value(index1267, GE_G[i, j]);
                        feature1.Store();
                    }
                }

                //保存法矢量线要素
                for (int i = 0; i < NodeCount; i++)
                {
                    for (int j = 0; j < TinCount - 1; j++)
                    {
                        if (NormalVector_F[i, j] == null)
                            continue;
                        //Focus
                        IFeature feature1 = NormalVectorLineFeatureClass_F.CreateFeature();
                        feature1.Shape = NormalVector_F[i, j] as IGeometry;
                        int index1 = feature1.Fields.FindField("Magnitude");
                        feature1.set_Value(index1, NLength_D);
                        int index167 = feature1.Fields.FindField("Probability");
                        feature1.set_Value(index167, 0.95 - (0.9 / 49) * j);
                        feature1.Store();
                    }
                }

                for (int i = 0; i < NodeCount; i++)
                {
                    for (int j = 0; j < TinCount - 1; j++)
                    {
                        if (NormalVector_D[i, j] == null)
                            continue;
                        //Disperse
                        IFeature feature10 = NormalVectorLineFeatureClass_D.CreateFeature();
                        feature10.Shape = NormalVector_D[i, j] as IGeometry;
                        int index2 = feature10.Fields.FindField("Magnitude");
                        feature10.set_Value(index2, NLength_D);
                        int index167 = feature10.Fields.FindField("Probability");
                        feature10.set_Value(index167, 0.95 - (0.9 / 49) * j);
                        feature10.Store();
                    }
                }

                //法矢量Focus箭头
                for (int i = 0; i < NodeCount; i++)
                {
                    for (int j = 0; j < TinCount - 1; j++)
                    {
                        IFeature feature11 = NormalVectorArrowFeatureClass_F.CreateFeature();
                        feature11.Shape = NormalArrow_F[i, j] as IMultiPatch;
                        feature11.Store();
                    }
                }

                //法矢量Disperse箭头
                for (int i = 0; i < NodeCount; i++)
                {
                    for (int j = 0; j < TinCount - 1; j++)
                    {
                        if (NormalArrow_D[i, j] == null)
                            continue;
                        IFeature feature11 = NormalVectorArrowFeatureClass_D.CreateFeature();
                        feature11.Shape = NormalArrow_D[i, j] as IMultiPatch;
                        feature11.Store();
                    }
                }

                //主方向           
                //创建主曲率矢量
                for (int i = 0; i < NodeCount; i++)
                {
                    for (int j = 0; j < TinCount - 1; j++)
                    {
                        if ((e1_F[i, j] == null) || (e2_F[i, j] == null))
                            continue;
                        //e1 Focus
                        IFeature feature6 = e1LineFeatureClass_F.CreateFeature();
                        feature6.Shape = e1_F[i, j] as IGeometry;
                        int index1 = feature6.Fields.FindField("k1");
                        feature6.set_Value(index1, k1_G[i, j]);
                        int index167 = feature6.Fields.FindField("Probability");
                        feature6.set_Value(index167, 0.95 - (0.9 / 49) * j);
                        feature6.Store();
                        //箭头
                        IFeature feature17 = e1ArrowFeatureClass_F.CreateFeature();
                        feature17.Shape = e1Arrow_F[i, j] as IMultiPatch;
                        if (k1_G[i, j] >= 0)
                        {
                            IsRegular = 1;
                        }
                        else
                        {
                            IsRegular = 0;
                        }
                        int index121 = feature17.Fields.FindField("IsRegular");
                        feature17.set_Value(index121, IsRegular);
                        feature17.Store();

                        //e1 Disperse
                        IFeature feature8 = e1LineFeatureClass_D.CreateFeature();
                        feature8.Shape = e1_D[i, j] as IGeometry;
                        int index11 = feature8.Fields.FindField("k1");
                        feature8.set_Value(index11, k1_G[i, j]);
                        int index267 = feature8.Fields.FindField("Probability");
                        feature8.set_Value(index267, 0.95 - (0.9 / 49) * j);
                        feature8.Store();
                        //箭头
                        IFeature feature71 = e1ArrowFeatureClass_D.CreateFeature();
                        feature71.Shape = e1Arrow_D[i, j] as IMultiPatch;
                        if (k1_G[i, j] >= 0)
                        {
                            IsRegular = 1;
                        }
                        else
                        {
                            IsRegular = 0;
                        }
                        int index01 = feature71.Fields.FindField("IsRegular");
                        feature71.set_Value(index01, IsRegular);
                        feature71.Store();

                        //e2 Focus
                        IFeature feature7 = e2LineFeatureClass_F.CreateFeature();
                        feature7.Shape = e2_F[i, j] as IGeometry;
                        int index2 = feature7.Fields.FindField("k2");
                        feature7.set_Value(index2, k2_G[i, j]);
                        int index1267 = feature7.Fields.FindField("Probability");
                        feature7.set_Value(index1267, 0.95 - (0.9 / 49) * j);
                        feature7.Store();
                        //箭头
                        IFeature feature23 = e2ArrowFeatureClass_F.CreateFeature();
                        feature23.Shape = e2Arrow_F[i, j] as IMultiPatch;
                        if (k2_G[i, j] >= 0)
                        {
                            IsRegular = 1;
                        }
                        else
                        {
                            IsRegular = 0;
                        }
                        int index21 = feature23.Fields.FindField("IsRegular");
                        feature23.set_Value(index21, IsRegular);
                        feature23.Store();

                        //e2 Disperse
                        IFeature feature2 = e2LineFeatureClass_D.CreateFeature();
                        feature2.Shape = e2_D[i, j] as IGeometry;
                        int index112 = feature2.Fields.FindField("k2");
                        feature2.set_Value(index112, k2_G[i, j]);
                        int index2267 = feature2.Fields.FindField("Probability");
                        feature2.set_Value(index2267, 0.95 - (0.9 / 49) * j);
                        feature2.Store();
                        //箭头
                        IFeature feature123 = e2ArrowFeatureClass_D.CreateFeature();
                        feature123.Shape = e2Arrow_D[i, j] as IMultiPatch;
                        if (k2_G[i, j] >= 0)
                        {
                            IsRegular = 1;
                        }
                        else
                        {
                            IsRegular = 0;
                        }
                        int index321 = feature123.Fields.FindField("IsRegular");
                        feature123.set_Value(index321, IsRegular);
                        feature123.Store();
                    }
                }

                //Normal首尾相连   
                IPoint[] SLstart = StartPoints;
                for (int p = 0; p < NodeCount; p++)
                {
                    for (int q = (TinCount - 1) / 2 - 1; q > -1; q--)
                    {
                        if (NormalVectors[p, q] == null)
                            continue;
                        IVector3D antivector = new Vector3DClass();
                        antivector.XComponent = -(NormalVectors[p, q].XComponent);
                        antivector.YComponent = -(NormalVectors[p, q].YComponent);
                        antivector.ZComponent = -(NormalVectors[p, q].ZComponent);
                        IRay ray = new RayClass();
                        ray.Origin = SLstart[p];
                        ray.Vector = antivector;
                        IPoint endPoint = new PointClass();
                        Material.MakeZAware(endPoint as IGeometry);
                        //设置矢量可视化长度
                        endPoint = ray.GetPointAtDistance(ConnectLength);                                                                                                   //cd 1
                        SLstart[p] = endPoint;
                    }
                }

                for (int p = 0; p < NodeCount; p++)
                {
                    //对每一个节点作流线
                    for (int q = 0; q < TinCount - 1; q++)
                    {
                        if (SLstart[p] == null)
                            continue;
                        if (NormalVectors[p, q] == null)
                            continue;
                        IRay ray = new RayClass();
                        ray.Origin = SLstart[p];
                        ray.Vector = NormalVectors[p, q];
                        IPoint endPoint = new PointClass();
                        Material.MakeZAware(endPoint as IGeometry);
                        //设置矢量可视化长度
                        endPoint = ray.GetPointAtDistance(ConnectLength);                                                                                                   //cd 1                                  
                        //添加两点构成向量的线端点
                        IPointCollection VectorPointCollection = new PolylineClass();
                        Material.MakeZAware(VectorPointCollection as IGeometry);
                        VectorPointCollection.AddPoint(SLstart[p], ref _missing, ref _missing);
                        VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                        //为向量末端添加箭头符号
                        IGeometry Arrow = Material.GetArrow();
                        double Inclination = NormalVectors[p, q].Inclination;
                        double degreesOfRotation = Math.Acos((NormalVectors[p, q].DotProduct(VerticalVector)) / ((NormalVectors[p, q].Magnitude) * (VerticalVector.Magnitude)));
                        ITransform3D transform3D = Arrow as ITransform3D;
                        //放缩
                        transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                        //转动
                        if (degreesOfRotation != 0)
                        {
                            double angleOfRotationInRadians = degreesOfRotation;
                            IVector3D axisOfRotationVector3D = new Vector3DClass();
                            axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, NormalVectors[p, q]);
                            transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                        }
                        //平移
                        if (endPoint.IsEmpty)
                            continue;
                        transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                        //创建要素
                        IFeature SLfeature2 = ConnectSLFeatureClass.CreateFeature();
                        SLfeature2.Shape = VectorPointCollection as IGeometry;
                        int index2267 = SLfeature2.Fields.FindField("Probability");
                        SLfeature2.set_Value(index2267, 0.95 - (0.9 / 49) * q);
                        SLfeature2.Store();
                        IFeature SLfeature3 = ConnectSLArrowFeatureClass.CreateFeature();
                        SLfeature3.Shape = Arrow as IMultiPatch;
                        SLfeature3.Store();
                        SLstart[p] = endPoint;
                    }
                }

                //********************************************************************************************
                //                                              写gocad的txt文件，分面存点
                //********************************************************************************************
                //曲率及其曲率熵
                StreamWriter MFile = new StreamWriter(MFilepath);
                StreamWriter MEFile = new StreamWriter(MEFilepath);
                StreamWriter GFile = new StreamWriter(GFilepath);
                StreamWriter GEFile = new StreamWriter(GEFilepath);

                for (int j = 0; j < TinCount; j++)
                {
                    string ZNodeClass = "ZNodeClass" + j;
                    IFields ZNodefields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPoint, ZNodeClass);
                    IFeatureClass ZNodeFeatureClass = CreateStandaloneFeatureClass(Workspace, ZNodeClass, ZNodefields);
                    //流
                    StreamWriter PointsFile = new StreamWriter(eachgocadpath + j + ".vs");
                    bool Head = true;
                    for (int i = 0; i < NodeCount; i++)
                    {
                        if ((Points_All[i, j] == null) || (Points_All[i, j].IsEmpty))
                        {
                            continue;
                        }

                        //所有点集
                        IFeature feature1 = ZNodeFeatureClass.CreateFeature();
                        feature1.Shape = Points_All[i, j];
                        int index11 = feature1.Fields.FindField("X");
                        feature1.set_Value(index11, Points_All[i, j].X);
                        int index12 = feature1.Fields.FindField("Y");
                        feature1.set_Value(index12, Points_All[i, j].Y);
                        int index13 = feature1.Fields.FindField("Z");
                        feature1.set_Value(index13, Points_All[i, j].Z);
                        int index3 = feature1.Fields.FindField("BigCurvature");
                        feature1.set_Value(index3, k1_G[i, j]);
                        int index4 = feature1.Fields.FindField("SmallCurvature");
                        feature1.set_Value(index4, k2_G[i, j]);
                        int index1 = feature1.Fields.FindField("GaussCurvature");
                        feature1.set_Value(index1, G_G[i, j]);
                        int index2 = feature1.Fields.FindField("MeanCurvature");
                        feature1.set_Value(index2, M_G[i, j]);
                        int index167 = feature1.Fields.FindField("Probability");
                        if (j == 50)
                        {
                            feature1.set_Value(index167, 0.5);
                        }
                        else
                        {
                            feature1.set_Value(index167, 0.95 - (0.9 / 49) * j);
                        }
                        int index267 = feature1.Fields.FindField("MEntropy");
                        feature1.set_Value(index267, ME_G[i, j]);
                        int index1267 = feature1.Fields.FindField("GEntropy");
                        feature1.set_Value(index1267, GE_G[i, j]);
                        feature1.Store();

                        //写入vs文件
                        string vsX = (Points_All[i, j].X).ToString();
                        string vsY = (Points_All[i, j].Y).ToString();
                        string vsZ = (Points_All[i, j].Z).ToString();
                        string vsG = (G_G[i, j]).ToString();
                        string vsM = (M_G[i, j]).ToString();
                        string vsME = (ME_G[i, j]).ToString();
                        string vsGE = (GE_G[i, j]).ToString();
                        string vsP;
                        if (j == 50)
                        {
                            vsP = (0.5).ToString();
                        }
                        else
                        {
                            vsP = (0.95 - (0.9 / 49) * j).ToString();
                        }
                        string vsOS = (GetOriginal.OSTD[i]).ToString();
                        if (Head == true)
                        {
                            PointsFile.WriteLine("GOCAD VSet" + "\r\n"
                                + "HEADER" + "\r\n"
                                + "{name: points" + j + "}" + "\r\n"
                                + "GOCAD_ORIGINAL_COORDINATE_SYSTEM\nNAME" + "\r\n"
                                + "PROJECTION Unknown" + "\r\n"
                                + "DATUM Unknown" + "\r\n"
                                + "AXIS_NAME X Y Z" + "\r\n"
                                + "AXIS_UNIT m m m" + "\r\n"
                                + "ZPOSITIVE Elevation" + "\r\n"
                                + "END_ORIGINAL_COORDINATE_SYSTEM" + "\r\n"
                                + "PROPERTIES G H HE GE P OS" + "\r\n"
                                + "PROP_LEGAL_RANGES **none**  **none** **none**  **none** **none**  **none** **none**  **none** **none**  **none** **none**  **none**" + "\r\n"
                                + "NO_DATA_VALUES -99999 -99999 -99999 -99999 -99999 -99999" + "\r\n"
                                + "PROPERTY_CLASSES g h he ge p os" + "\r\n"
                                + "PROPERTY_KINDS \"Real Number\" \"Real Number\" \"Real Number\" \"Real Number\" \"Real Number\" \"Real Number\"" + "\r\n"
                                + "PROPERTY_SUBCLASSES QUANTITY Float QUANTITY Float QUANTITY Float QUANTITY Float QUANTITY Float QUANTITY Float" + "\r\n"
                                + "ESIZES 1  1  1  1  1  1" + "\r\n"
                                + "UNITS unitless unitless unitless unitless unitless unitless " + "\r\n"
                                + "PROPERTY_CLASS_HEADER X {" + "\r\n" + "kind: X" + "\r\n" + "unit: m" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER Y {" + "\r\n" + "kind: Y" + "\r\n" + "unit: m" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER Z {" + "\r\n" + "kind: Depth\nunit: m" + "\r\n" + "is_z: on" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER g {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER h {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER he {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER ge {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER p {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER os {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                                );
                            Head = false;
                        }
                        PointsFile.Write("PVRTX " + i + " ");
                        PointsFile.Write(vsX + " ");
                        PointsFile.Write(vsY + " ");
                        PointsFile.Write(vsZ + " ");
                        PointsFile.Write(vsG + " ");
                        PointsFile.Write(vsM + " ");
                        PointsFile.Write(vsME + " ");
                        PointsFile.Write(vsGE + " ");
                        PointsFile.Write(vsP + " ");
                        PointsFile.WriteLine(vsOS + " ");
                        //写入H，HE，G，GE
                        MFile.Write(vsM + " ");
                        GFile.Write(vsG + " ");
                        MEFile.Write(vsME + " ");
                        GEFile.Write(vsGE + " ");
                    }
                    MFile.WriteLine("");
                    MEFile.WriteLine("");
                    GFile.WriteLine("");
                    GEFile.WriteLine("");
                    PointsFile.Write("END");
                    PointsFile.Close();
            #endregion Creat Global Feature
                }
                MFile.Close();
                MEFile.Close();
                GFile.Close();
                GEFile.Close();

                //输出法矢量的XYZ坐标值
                StreamWriter NFile = new StreamWriter(NVXYZpath);
                for (int i = 0; i < TinCount; i++)
                {
                    for (int j = 0; j < NodeCount - 4; j++)
                    {
                        NFile.Write(NormalVectors[j, i].XComponent + " " + NormalVectors[j, i].YComponent + " " + NormalVectors[j, i].ZComponent + ",");
                    }
                    NFile.WriteLine("");
                }
                NFile.Close();

                //输出e1矢量的XYZ坐标值
                StreamWriter E1File = new StreamWriter(E1XYZpath);
                for (int i = 0; i < TinCount; i++)
                {
                    for (int j = 0; j < NodeCount - 4; j++)
                    {
                        E1File.Write(e1s[j, i].XComponent + " " + e1s[j, i].YComponent + " " + e1s[j, i].ZComponent + " ");
                    }
                    E1File.WriteLine("");
                }
                E1File.Close();

                //输出e2矢量的XYZ坐标值
                StreamWriter E2File = new StreamWriter(E2XYZpath);
                for (int i = 0; i < TinCount; i++)
                {
                    for (int j = 0; j < NodeCount - 4; j++)
                    {
                        E2File.Write(e2s[j, i].XComponent + " " + e2s[j, i].YComponent + " " + e2s[j, i].ZComponent + " ");
                    }
                    E2File.WriteLine("");
                }
                E2File.Close();
            }
        }

        #region 排序函数
        public int NodeLocation(ITinNode NodeL)
        {
            int v = (int)((NodeL.X - (XMin - (RasterSize / 2))) / RasterSize) + RowCount * (int)((NodeL.Y - (YMin - (RasterSize / 2))) / (RasterSize));
            return v;
        }

        public int PointLocation(IPoint PointL)
        {
            int v = (int)((PointL.X - (XMin - (RasterSize / 2))) / RasterSize) + RowCount * (int)((PointL.Y - (YMin - (RasterSize / 2))) / (RasterSize));
            return v;
        }
        #endregion
    }
}
        #endregion CreateFeature






