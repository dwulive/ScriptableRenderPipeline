using System.Reflection;

namespace UnityEngine.MaterialGraph
{
    [Title("Math/Round/Ceiling")]
    public class CeilingNode : CodeFunctionNode
    {
        public CeilingNode()
        {
            name = "Ceiling";
        }

        protected override MethodInfo GetFunctionToConvert()
        {
            return GetType().GetMethod("Unity_Ceiling", BindingFlags.Static | BindingFlags.NonPublic);
        }

        static string Unity_Ceiling(
            [Slot(0, Binding.None)] DynamicDimensionVector In,
            [Slot(1, Binding.None)] out DynamicDimensionVector Out)
        {
            return
                @"
{
    Out = ceil(In);
}
";
        }
    }
}
