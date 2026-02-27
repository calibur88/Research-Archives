import requests

def download_ecdata_chunk(conductor_start, conductor_end):
    """
    下载指定导子范围的 allbsd 数据块
    
    例如：导子 0-9999 对应文件 allbsd.00000-09999（无扩展名）
    """
    # 格式化文件名 (例如 00000-09999)，注意：没有 .txt 扩展名！
    start_str = f"{conductor_start:05d}"
    end_str = f"{conductor_end:05d}"
    filename = f"allbsd.{start_str}-{end_str}"  # 修正：去掉 .txt
    
    # ecdata 仓库的 raw 文件 URL
    base_url = "https://raw.githubusercontent.com/JohnCremona/ecdata/master/allbsd"
    file_url = f"{base_url}/{filename}"
    
    print(f"正在下载: {file_url}")
    response = requests.get(file_url)
    
    if response.status_code == 200:
        # 保存到本地（可选择是否添加 .txt 扩展名便于本地识别）
        local_filename = f"{filename}.txt"
        with open(local_filename, 'w', encoding='utf-8') as f:
            f.write(response.text)
        print(f"成功保存: {local_filename}")
        return response.text
    else:
        print(f"下载失败，状态码: {response.status_code}")
        print(f"请检查导子范围 {conductor_start}-{conductor_end} 是否存在对应文件")
        return None

# 使用示例：下载导子 0-9999 的数据
data = download_ecdata_chunk(0, 9999)
