using UnityEngine.Rendering.HDPipeline;

namespace UnityEditor.Rendering.HDPipeline
{
    // Editor only functions for HDAdditoonalLightData User API
    public static class HDAdditionalLightDataEditorExtension
    {
        /// <summary>
        /// Toggle the usage of color temperature.
        /// </summary>
        /// <param name="hdLight"></param>
        /// <param name="enable"></param>
        public static void EnableColorTemperature(this HDAdditionalLightData hdLight, bool enable)
        {
            hdLight.useColorTemperature = enable;
        }
    }
}
