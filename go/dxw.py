# -*- coding: utf-8 -*-

from django.shortcuts import render
from django.views.decorators import csrf
from pc import daoxian

info1 = """
 0 -- 1	65.893
 1 -- 2	103.656
 2 -- 3	96.446
 3 -- 4	82.557
 4 -- 5	已知值
 5 -- 6	71.194
 6 -- 0	102.801
 0 -- 7	63.84
 7 -- 8	70.315
 8 -- 3	47.726
 """

info2 = """
∠0 1 2	231°50′15.1″
∠1 0 6	134°25′21.9″
∠1 0 7	62°25′37.2″
∠1 2 3	251°31′27.3″
∠2 3 4	204°32′33.3″
∠2 3 8	288°18′12.4″
∠3 4 5	246°45′15.2″
∠4 3 8	83°45′53.7″
∠4 5 6	244°09′57.8″
∠5 6 0	215°35′53.8″
∠6 0 7	288°00′32″
∠0 7 8	192°59′25.3″
∠7 8 3	156°15′13.5″
"""

info3 = '''已知点：4 x=97.478(m), y=155.682(m)
已知点：5 x=184.833(m), y=90.147(m)
测边误差：1/2000
测角误差：12"'''

weight_ratio = '1'


# 接收POST请求数据
def cal_dxw(request):
    ctx = {
        'info1': info1, 'info2': info2, 'info3': info3, 'weight_ratio': weight_ratio
    }
    if request.POST:
        ctx['params'] = 'producing'
        try:
            ctx['info1'], ctx['info2'], ctx['info3'], ctx['weight_ratio'] = request.POST['info1'], request.POST['info2'], request.POST['info3'], request.POST['weight_ratio']
            dxw = daoxian.Daoxian(ctx['info1'], ctx['info2'], ctx['info3'], ctx['weight_ratio'])
            dxw.check_info()
            dxw.init_params()
            dxw.init_weight()
            dxw.cal()
            ctx['adjust_result'] = '<h2>平差结果</h2>'
            ctx['params'] = '<p>params</p>' + dxw.html_of_params()
            ctx['edges'] = '<p>edges observation</p>' + dxw.html_of_edges()
            ctx['corners'] = '<p>corners observation</p>' + dxw.html_of_corners()
            ctx['edge_sigma_max'] = '<p>平差后最弱边中误差（mm）</p>' + dxw.error['edge_sigma_max'].__str__()
            ctx['corner_sigma_max'] = '<p>平差后最弱角中误差（″）</p>' + dxw.error['corner_sigma_max'].__str__()
            ctx['point_sigma_max'] = '<p>平差后最弱点位中误差（mm）</p>' + dxw.error['point_sigma_max'].__str__()
            ctx['sigma_final'] = '<p>平差后单位权中误差（″）</p>' + dxw.error['sigma_final'].__str__()

            ctx['intermediates'] = '<h2>中间变量</h2>'
            ctx['b_coefficient'] = '<p>系数矩阵 B. 其中角度单位为″，长度单位为mm</p>' + dxw.intermediates['b_coefficient']
            ctx['l_observation_residual'] = '<p>辅助量 l. 其中角度单位为″，长度单位为mm</p>' + dxw.intermediates['l_observation_residual']
            ctx['v_residual'] = '<p>观测值改正数 v. 其中角度单位为″，长度单位为mm</p>' + dxw.intermediates['v_residual']
            ctx['P_weight_matrix'] = '<p>权重矩阵 P</p>' + dxw.intermediates['P_weight_matrix']
        except:
            ctx['adjust_result'] = '<p>error input</p>'
    return render(request, "dxw_solve.html", ctx)