using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using Grasshopper.Kernel;
using Rhino.Geometry;
using ClosedXML.Excel;
using MathNet.Numerics;

namespace SolidFEM.Components
{
    public class WriteListToColumn : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public WriteListToColumn()
          : base("ListToColumn", "lst2col",
              "Write a column to a list in an Excel file.",
              "SmartMesh", "Data")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Path", "p", "Path to Excel file", GH_ParamAccess.item); // 0
            pManager.AddTextParameter("Sheet", "s", "Sheet to write to", GH_ParamAccess.item); // 1
            pManager.AddTextParameter("Column", "col", "Column to write to", GH_ParamAccess.item); // 2
            pManager.AddNumberParameter("DataList", "list", "List of data to write", GH_ParamAccess.list); // 3
            pManager.AddBooleanParameter("Write", "w", "Set to true to write to file", GH_ParamAccess.item, false); // 4

            pManager[1].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // empty output
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // -- input --
            string path = ""; // 0
            string sheet = ""; // 1
            string col = ""; // 2
            List<double> dataLst = new List<double>(); // 3
            bool run = false; // 4

            if (!DA.GetData(0, ref path))return; // 0
            DA.GetData(1, ref sheet); // 1
            DA.GetData(2, ref col); // 2
            if (!DA.GetDataList(3, dataLst))return ; // 3
            DA.GetData(4, ref run); // 4
            // -- solve --

            if (!run)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Set Write to 'true' to run component");
                return; 
            }

            // try to open the fil
            XLWorkbook wb = new XLWorkbook();
            try
            {
                wb = new XLWorkbook(path);
            }
            catch (Exception e)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, e.ToString());
            }

            IXLWorksheet ws;
            if (!wb.TryGetWorksheet(sheet, out ws))
            {
                ws = wb.Worksheet(0); // if the worksheet does not exist, use the standard one. 
            }

            // cast col to int or string. Currently not working... Need some control here I guess...


            if (!Int32.TryParse(col, out int colInt))
            {
                // it is a string
                for (int i = 0; i < dataLst.Count; i++)
                {
                    ws.Column(col).Cell(i + 2).Value = dataLst[i]; // add +2 to row to avoid overwriting the header (Excel index starting from 1).
                }
            }
            else
            {
                // the column is an integer
                // write data to specified column. 
                for (int i = 0; i < dataLst.Count; i++)
                {
                    ws.Column(colInt).Cell(i + 2).Value = dataLst[i]; // add +1 to row to avoid overwriting the header.
                }
            }

            wb.SaveAs(path);
            /*
            string name = "Excel Fra Component.xlsx";
            //var path = Directory.GetCurrentDirectory();
            var path2 = "C:\\Users\\sverremh\\OneDrive - NTNU\\Spring 2022\\Special Issue Solid FEM";
            IXLWorkbook wb = new XLWorkbook();
            IXLWorksheet ws = wb.Worksheets.Add("Sample Sheet");

            ws.Cell(1, 1).Value = "Hello World With new framework !";
            ws.Cell(2, 1).Value = "Hello From Grasshopper !";

            wb.SaveAs(path2 + "\\" + name);
            // wb.SaveAs(path + "\\" + name);
            */

            // -- output
            AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Column of info added to file");
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("2757e60a-e084-4594-9200-7ee395ea8f81"); }
        }
    }
}