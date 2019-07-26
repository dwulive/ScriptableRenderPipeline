using System;
using u;
using UnityEngine;

namespace u
{
    public struct I
    {
        public Vector3 camC;
        public Matrix4x4 projMatrix0;
        public Matrix4x4 projMatrix1;
        public float nopt;

        public static float nearZ = 1.0f;
        public static float farZ = 128.0f;

        public static float nearZ0 = 0.5f;
        public static float farZ0 = 1024;

        public static float nearZ1 = 1.0f;
        public static float farZ1 = 128.0f;

        public static float yGain = 1.0f;
        public static float bGain = 1.0f;
        public static float wGain = 1;
        public static float aGain = 1;
        public static float xGain = 1;
        public static float pGain = -1;
        public static bool swapW = true;
        public static bool pOrder = false;
        public static float minH;
        public static float camYGain = -1;

        public static I i;
    }

}
namespace UnityEngine.Rendering.LWRP
{
    public struct ShadowSliceData
    {
        public Matrix4x4 viewMatrix;
        public Matrix4x4 projectionMatrix;
        public Matrix4x4 shadowTransform;
        public int offsetX;
        public int offsetY;
        public int resolution;

        public void Clear()
        {
            viewMatrix = Matrix4x4.identity;
            projectionMatrix = Matrix4x4.identity;
            shadowTransform = Matrix4x4.identity;
            offsetX = offsetY = 0;
            resolution = 1024;
        }
    }

    public static class ShadowUtils
    {
        public static  RenderTextureFormat m_ShadowmapFormat;
        public static  bool m_ForceShadowPointSampling;

        static ShadowUtils()
        {
            m_ShadowmapFormat = RenderingUtils.SupportsRenderTextureFormat(RenderTextureFormat.Shadowmap) && (SystemInfo.graphicsDeviceType != GraphicsDeviceType.OpenGLES2)
                ? RenderTextureFormat.Shadowmap
                : RenderTextureFormat.Depth;
            m_ForceShadowPointSampling = SystemInfo.graphicsDeviceType == GraphicsDeviceType.Metal &&
                GraphicsSettings.HasShaderDefine(Graphics.activeTier, BuiltinShaderDefine.UNITY_METAL_SHADOWS_USE_POINT_FILTERING);
        }

        public static bool ExtractDirectionalLightMatrix(float LoV,Vector3 camZ,Vector3 camC,ref CullingResults cullResults, ref ShadowData shadowData, int shadowLightIndex, int cascadeIndex, int shadowmapWidth, int shadowmapHeight, int shadowResolution, float shadowNearPlane, out Vector4 cascadeSplitDistance, out ShadowSliceData shadowSliceData, out Matrix4x4 viewMatrix, out Matrix4x4 projMatrix)
        {
            ShadowSplitData splitData;
            bool success = cullResults.ComputeDirectionalShadowMatricesAndCullingPrimitives(shadowLightIndex,
                cascadeIndex, shadowData.mainLightShadowCascadesCount, shadowData.mainLightShadowCascadesSplit, shadowResolution, shadowNearPlane, out viewMatrix, out projMatrix,
                out splitData);
            var temp = viewMatrix.inverse;
            var y0 = (Vector3)temp.GetColumn(1);
            var dot0 = Vector3.Dot(y0, camZ);
            var z0 = (Vector3)temp.GetColumn(2);
            var camZ0 =camZ- Vector3.Dot(camZ, z0) * z0;
            camZ0.Normalize();
            var cr = Vector3.Cross(y0, camZ0);
            var si = Vector3.Dot(cr, z0);
            var cs = Vector3.Dot(y0, camZ0);
            var theta = Mathf.Atan2(si, cs);
            var rot = Matrix4x4.Rotate(Quaternion.Euler(0, 0, -theta * Mathf.Rad2Deg));
            temp = temp * rot;
            var y1 = (Vector3)temp.GetColumn(1);
            var dot1 = Vector3.Dot(y1, camZ);
            var dot2 = Vector3.Dot(y1, camZ0);

            viewMatrix = temp.inverse;

            applyLISPSM(LoV, camC, ref viewMatrix, ref projMatrix);
            cascadeSplitDistance = splitData.cullingSphere;
            shadowSliceData.offsetX = (cascadeIndex % 2) * shadowResolution;
            shadowSliceData.offsetY = (cascadeIndex / 2) * shadowResolution;
            shadowSliceData.resolution = shadowResolution;
            shadowSliceData.viewMatrix = viewMatrix;
            shadowSliceData.projectionMatrix = projMatrix;
            shadowSliceData.shadowTransform = GetShadowTransform(projMatrix, viewMatrix);

            // If we have shadow cascades baked into the atlas we bake cascade transform
            // in each shadow matrix to save shader ALU and L/S
            if (shadowData.mainLightShadowCascadesCount > 1)
                ApplySliceTransform(ref shadowSliceData, shadowmapWidth, shadowmapHeight);

            return success;
        }

        public static bool ExtractSpotLightMatrix(ref CullingResults cullResults, ref ShadowData shadowData, int shadowLightIndex, out Matrix4x4 shadowMatrix, out Matrix4x4 viewMatrix, out Matrix4x4 projMatrix)
        {
            ShadowSplitData splitData;
            bool success = cullResults.ComputeSpotShadowMatricesAndCullingPrimitives(shadowLightIndex, out viewMatrix, out projMatrix, out splitData);
            shadowMatrix = GetShadowTransform(projMatrix, viewMatrix);
            return success;
        }

        public static void RenderShadowSlice(CommandBuffer cmd, ref ScriptableRenderContext context,
            ref ShadowSliceData shadowSliceData, ref ShadowDrawingSettings settings,
            Matrix4x4 proj, Matrix4x4 view)
        {
            cmd.SetViewport(new Rect(shadowSliceData.offsetX, shadowSliceData.offsetY, shadowSliceData.resolution, shadowSliceData.resolution));
            cmd.EnableScissorRect(new Rect(shadowSliceData.offsetX + 4, shadowSliceData.offsetY + 4, shadowSliceData.resolution - 8, shadowSliceData.resolution - 8));

            cmd.SetViewProjectionMatrices(view, proj);
            context.ExecuteCommandBuffer(cmd);
            cmd.Clear();
            context.DrawShadows(ref settings);
            cmd.DisableScissorRect();
            context.ExecuteCommandBuffer(cmd);
            cmd.Clear();
        }

        public static int GetMaxTileResolutionInAtlas(int atlasWidth, int atlasHeight, int tileCount)
        {
            int resolution = Mathf.Min(atlasWidth, atlasHeight);
            int currentTileCount = atlasWidth / resolution * atlasHeight / resolution;
            while (currentTileCount < tileCount)
            {
                resolution = resolution >> 1;
                currentTileCount = atlasWidth / resolution * atlasHeight / resolution;
            }
            return resolution;
        }

        public static void ApplySliceTransform(ref ShadowSliceData shadowSliceData, int atlasWidth, int atlasHeight)
        {
            Matrix4x4 sliceTransform = Matrix4x4.identity;
            float oneOverAtlasWidth = 1.0f / atlasWidth;
            float oneOverAtlasHeight = 1.0f / atlasHeight;
            sliceTransform.m00 = shadowSliceData.resolution * oneOverAtlasWidth;
            sliceTransform.m11 = shadowSliceData.resolution * oneOverAtlasHeight;
            sliceTransform.m03 = shadowSliceData.offsetX * oneOverAtlasWidth;
            sliceTransform.m13 = shadowSliceData.offsetY * oneOverAtlasHeight;

            // Apply shadow slice scale and offset
            shadowSliceData.shadowTransform = sliceTransform * shadowSliceData.shadowTransform;
        }

        public static Vector4 GetShadowBias(ref VisibleLight shadowLight, int shadowLightIndex, ref ShadowData shadowData, Matrix4x4 lightProjectionMatrix, float shadowResolution)
        {
            if (shadowLightIndex < 0 || shadowLightIndex >= shadowData.bias.Count)
            {
                Debug.LogWarning(string.Format("{0} is not a valid light index.", shadowLightIndex));
                return Vector4.zero;
            }

            float frustumSize;
            if (shadowLight.lightType == LightType.Directional)
            {
                // Frustum size is guaranteed to be a cube as we wrap shadow frustum around a sphere
                frustumSize = 2.0f / lightProjectionMatrix.m00;
            }
            else if (shadowLight.lightType == LightType.Spot)
            {
                // For perspective projections, shadow texel size varies with depth
                // It will only work well if done in receiver side in the pixel shader. Currently LWRP
                // do bias on caster side in vertex shader. When we add shader quality tiers we can properly
                // handle this. For now, as a poor approximation we do a constant bias and compute the size of
                // the frustum as if it was orthogonal considering the size at mid point between near and far planes.
                // Depending on how big the light range is, it will be good enough with some tweaks in bias
                frustumSize = Mathf.Tan(shadowLight.spotAngle * 0.5f * Mathf.Deg2Rad) * shadowLight.range;
            }
            else
            {
                Debug.LogWarning("Only spot and directional shadow casters are supported in lightweight pipeline");
                frustumSize = 0.0f;
            }

            // depth and normal bias scale is in shadowmap texel size in world space
            float texelSize = frustumSize / shadowResolution;
            float depthBias = -shadowData.bias[shadowLightIndex].x * texelSize;
            float normalBias = -shadowData.bias[shadowLightIndex].y * texelSize;
            
            if (shadowData.supportsSoftShadows)
            {
                // TODO: depth and normal bias assume sample is no more than 1 texel away from shadowmap
                // This is not true with PCF. Ideally we need to do either
                // cone base bias (based on distance to center sample)
                // or receiver place bias based on derivatives.
                // For now we scale it by the PCF kernel size (5x5)
                const float kernelRadius = 2.5f;
                depthBias *= kernelRadius;
                normalBias *= kernelRadius;
            }

            return new Vector4(depthBias, normalBias, 0.0f, 0.0f);
        }

        public static void SetupShadowCasterConstantBuffer(CommandBuffer cmd, ref VisibleLight shadowLight, Vector4 shadowBias)
        {
            Vector3 lightDirection = -shadowLight.localToWorldMatrix.GetColumn(2);
            cmd.SetGlobalVector("_ShadowBias", shadowBias);
            cmd.SetGlobalVector("_LightDirection", new Vector4(lightDirection.x, lightDirection.y, lightDirection.z, 0.0f));
        }

        public static RenderTexture GetTemporaryShadowTexture(int width, int height, int bits)
        {
            var shadowTexture = RenderTexture.GetTemporary(width, height, bits, m_ShadowmapFormat);
            shadowTexture.filterMode = m_ForceShadowPointSampling ? FilterMode.Point : FilterMode.Bilinear;
            shadowTexture.wrapMode = TextureWrapMode.Clamp;

            return shadowTexture;
        }

        // This construct a frustum (similar to glFrustum or frustum), except
        // it looks towards the +y axis, and assumes -1,1 for the left/right and bottom/top planes.
 static  Matrix4x4 warpFrustum(float n, float f)
 {
         float d = 1 / (f - n);
         float A = (f + n) * d;
         float B = -2 * n * f * d;
         return new Matrix4x4 (
             new Vector4(n, 0, 0, 0),
             new Vector4(0, A* I.aGain, 0, I.swapW ? I.wGain : I.yGain * I.bGain *B),
             new Vector4(0, 0, n, 0),
             new Vector4(0, I.swapW ? I.yGain * I.bGain * B : I.wGain, 0, 0));
}

        
        static void applyLISPSM(float LoV, Vector3 camWC, ref Matrix4x4 viewMatrix, ref Matrix4x4 projMatrix)
{

 //    float LoV = dot(camera.getForwardVector(), dir);
     float sinLV = Mathf.Sqrt(1.0f - LoV * LoV);

            // Virtual near plane -- the default is 1m, can be changed by the user.
            // The virtual near plane prevents too much resolution to be wasted in the area near the eye
            // where shadows might not be visible (e.g. a character standing won't see shadows at her feet).
            //     float dzn = Mathf.Max(0.0f, params.options.shadowNearHint - camera.zn);
            //     float dzf = Mathf.Max(0.0f, camera.zf - params.options.shadowFarHint);
            float dzn = I.nearZ - I.nearZ0;// Mathf.Max(0.0f, params.options.shadowNearHint - camera.zn);
            float dzf = I.farZ0 - I.farZ;// Mathf.Max(0.0f, camera.zf - params.options.shadowFarHint);

            // near/far plane's distance from the eye in view space of the shadow receiver volume.
            //      Vector2 znf = -computeNearFar(camera.view, wsShadowReceiversVolume.data(), vertexCount);
            float zn = I.nearZ0;// Mathf.Max(camera.zn, znf[0]); // near plane distance from the eye
            float zf = I.farZ0;// Mathf.Min(camera.zf, znf[1]); // far plane distance from the eye

            Vector3 lsCameraPosition = viewMatrix * camWC;
            I.i.camC = lsCameraPosition;
            // compute n and f, the near and far planes coordinates of Wp (warp space).
            // It's found by looking down the Y axis in light space (i.e. -Z axis of Wp,
            // i.e. the axis orthogonal to the light direction) and taking the Min/Max
            // of the shadow receivers volume.
            // Note: znear/zfar encoded in Mp has no influence here (b/c we're interested only by the y axis)
            // Vector2 nf = computeNearFarOfWarpSpace(LMpMv, wsShadowReceiversVolume.data(), vertexCount);
            float n = I.camYGain *lsCameraPosition.y + I.nearZ1;// nf[0];              // near plane coordinate of Mp (light space)
            float f = n+ I.farZ1;//conservative estimate, we might be able to get away with less // nf[1];              // far plane coordinate of Mp (light space)
            float d = Mathf.Abs(f - n);    // Wp's depth-range d (abs necessary because we're dealing with z-coordinates, not distances)

    // The simplification below is correct only for directional lights
     float z0 = zn;                // for directional lights, z0 = zn
     float z1 = z0 + d * sinLV;    // btw, note that z1 doesn't depend on zf


          //  Matrix4x4 W = Matrix4x4.identity;
            
    // see nopt1 below for an explanation about this test
    if (3.0f * (dzn / (zf - zn)) < 2.0f)
    {
        // nopt is the optimal near plane distance of Wp (i.e. distance from P).

        // virtual near and far planes
         float vz0 = Mathf.Max(0.0f, Mathf.Max(zn + dzn, z0));
         float vz1 = Mathf.Max(0.0f, Mathf.Min(zf - dzf, z1));

        // in the general case, nopt is computed as:
         float nopt0 = (1.0f / sinLV) * (z0 + Mathf.Sqrt(vz0 * vz1) /* zero */ );

        // however, if dzn becomes too large, the Max error doesn't happen in the depth range,
        // and the equation below should be used instead. If dzn reaches 2/3 of the depth range
        // zf-zn, nopt becomes infinite, and we must revert to an ortho projection.
         float nopt1 = dzn / (2.0f - 3.0f * (dzn / (zf - zn)));

        // We simply use the Max of the two expressions
          var nopt = Mathf.Max(nopt0, nopt1);
          I.i.nopt = nopt;

         Vector3 p = new Vector3(
                // Another option here is to use lsShadowReceiversCenter.x, which skews less the
                // x axis. Doesn't seem to make a big difference in the end.
                lsCameraPosition.x* I.xGain,
                n - nopt,
                // note: various papers suggest to use the shadow receiver's center z coordinate in light
                // space, i.e. to center "vertically" on the shadow receiver volume.
                // e.g. (LMpMv * wsShadowReceiversVolume.center()).z
                // However, simply using 0, guarantees to be centered on the light frustum, which itself
                // is built from the shadow receiver and/or casters bounds.
                0
        );

         Matrix4x4 Wv = Matrix4x4.Translate(I.pGain *p);
         viewMatrix = Wv * viewMatrix;
                var Wp = warpFrustum(nopt, nopt + d);
                I.i.projMatrix0 = projMatrix;
                I.i.projMatrix1 = Wp;
                projMatrix = I.pOrder ? projMatrix*Wp : Wp * projMatrix;
               // W = Wp;
 //               W = Wp * Wv;
    }
            else
            {
                Debug.Log("Invalid D");

            }

}

static Matrix4x4 GetShadowTransform(Matrix4x4 proj, Matrix4x4 view)
        {
            // Currently CullResults ComputeDirectionalShadowMatricesAndCullingPrimitives doesn't
            // apply z reversal to projection matrix. We need to do it manually here.
            if (SystemInfo.usesReversedZBuffer)
            {
                proj.m20 = -proj.m20;
                proj.m21 = -proj.m21;
                proj.m22 = -proj.m22;
                proj.m23 = -proj.m23;
            }

            Matrix4x4 worldToShadow = proj * view;

            var textureScaleAndBias = Matrix4x4.identity;
            textureScaleAndBias.m00 = 0.5f;
            textureScaleAndBias.m11 = 0.5f;
            textureScaleAndBias.m22 = 0.5f;
            textureScaleAndBias.m03 = 0.5f;
            textureScaleAndBias.m23 = 0.5f;
            textureScaleAndBias.m13 = 0.5f;

            // Apply texture scale and offset to save a MAD in shader.
            return textureScaleAndBias * worldToShadow;
        }
    }
}
