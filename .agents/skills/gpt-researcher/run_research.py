#!/usr/bin/env python3
"""
GPT Researcher 包装器脚本 — Conda 环境版本

用法:
    python .agents/skills/gpt-researcher/run_research.py "你的研究问题"
    python .agents/skills/gpt-researcher/run_research.py --report-type deep "深度研究问题"
    python .agents/skills/gpt-researcher/run_research.py --output my_report.md "问题"
    python .agents/skills/gpt-researcher/run_research.py --check-keys
"""

import argparse
import asyncio
import os
import re
import sys
from datetime import datetime
from pathlib import Path

# ── 项目根目录 ──────────────────────────────────────────────────
_PROJECT_ROOT = Path(__file__).resolve().parents[3]  # .agents/skills/gpt-researcher/ → root

# ── 将 conda site-packages 加入 sys.path ─────────────────────────
import site
_VENV_SITE_PACKAGES = site.getsitepackages()[0]
if _VENV_SITE_PACKAGES not in sys.path:
    sys.path.insert(0, _VENV_SITE_PACKAGES)

# ── 加载 .env ───────────────────────────────────────────────────
try:
    from dotenv import load_dotenv
    _env_path = _PROJECT_ROOT / ".env"
    if _env_path.exists():
        load_dotenv(_env_path)
        print(f"[dotenv] Loaded: {_env_path}")
    else:
        print(f"[dotenv] WARNING: .env not found at {_env_path}")
except ImportError:
    print("[dotenv] WARNING: python-dotenv not installed; skipping .env loading")


def check_api_keys() -> bool:
    """验证所有必需的 API Key 是否已设置"""
    required = {
        "DEEPSEEK_API_KEY": os.getenv("DEEPSEEK_API_KEY"),
        "TAVILY_API_KEY": os.getenv("TAVILY_API_KEY"),
    }
    all_set = True
    for name, value in required.items():
        if value:
            print(f"  ✅ {name}: {'*' * 20}{value[-4:] if len(value) > 4 else '****'}")
        else:
            print(f"  ❌ {name}: MISSING")
            all_set = False

    # 可选
    embedding = os.getenv("EMBEDDING", "")
    if embedding:
        print(f"  ✅ EMBEDDING: {embedding}")
    else:
        print(f"  ⚠️  EMBEDDING: not set (will use default)")

    return all_set


def sanitize_filename(query: str, max_len: int = 50) -> str:
    """将查询转为安全的文件名片段"""
    # 取前 max_len 个非空白字符，组成 slug
    slug = re.sub(r'[^\w一-鿿-]', '_', query.strip())
    slug = re.sub(r'_+', '_', slug).strip('_')
    if len(slug) > max_len:
        slug = slug[:max_len].rstrip('_')
    return slug


def default_output_path(query: str, report_type: str) -> Path:
    """生成默认输出路径: references/{YYMMDD}_{query摘要}.md"""
    reports_dir = _PROJECT_ROOT / "references"
    reports_dir.mkdir(parents=True, exist_ok=True)
    date_prefix = datetime.now().strftime("%y%m%d")
    slug = sanitize_filename(query)
    # 如果 report_type 不是默认的 research_report，加入文件名
    if report_type != "research_report":
        slug = f"{report_type}_{slug}" if slug else report_type
    filename = f"{date_prefix}_{slug}.md" if slug else f"{date_prefix}_report.md"
    return reports_dir / filename


async def run_research(
    query: str,
    report_type: str = "research_report",
    output_path: Path | None = None,
    verbose: bool = True,
) -> str:
    """执行 GPT Researcher 研究并保存报告"""
    from gpt_researcher import GPTResearcher

    if output_path is None:
        output_path = default_output_path(query, report_type)

    if verbose:
        print(f"\n{'='*60}")
        print(f"🔬 GPT Researcher")
        print(f"{'='*60}")
        print(f"  Query:       {query}")
        print(f"  Report Type: {report_type}")
        print(f"  Output:      {output_path}")
        print(f"  LLM:         {os.getenv('FAST_LLM', 'deepseek:deepseek-v4-pro')}")
        print(f"{'='*60}\n")

    researcher = GPTResearcher(
        query=query,
        report_type=report_type,
        report_source="web",
        verbose=verbose,
    )

    # 执行研究
    print("⏳ Conducting research...")
    await researcher.conduct_research()

    # 撰写报告
    print("📝 Writing report...")
    report = await researcher.write_report()

    # 保存
    output_path.write_text(report, encoding="utf-8")
    print(f"\n✅ Report saved to: {output_path}")
    print(f"   Length: {len(report)} chars / ~{len(report.split())} words")

    return report


def main():
    parser = argparse.ArgumentParser(
        description="GPT Researcher — 包装器脚本 (Conda 环境)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s "What are the latest AI developments?"
  %(prog)s --report-type deep "Quantum computing advances"
  %(prog)s --output my_report.md "CRISPR gene editing"
  %(prog)s --check-keys
  %(prog)s --verbose "Research question"
        """,
    )
    parser.add_argument(
        "query",
        nargs="?",
        help="研究问题 / 搜索查询",
    )
    parser.add_argument(
        "--report-type",
        default="research_report",
        choices=[
            "research_report",
            "detailed_report",
            "deep",
            "outline_report",
            "resource_report",
        ],
        help="报告类型 (default: research_report)",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="自定义输出文件名 (自动保存到 references/)",
    )
    parser.add_argument(
        "--check-keys",
        action="store_true",
        help="验证 API Key 配置是否正确",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=True,
        help="详细输出 (default: True)",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="静默模式",
    )

    args = parser.parse_args()

    # ── check-keys 模式 ──
    if args.check_keys:
        print("🔑 Checking API Keys...")
        ok = check_api_keys()
        if ok:
            print("\n✅ All required API keys are set.")
        else:
            print("\n❌ Some API keys are missing. Please check your .env file.")
        sys.exit(0 if ok else 1)

    # ── 需要 query ──
    if not args.query:
        parser.print_help()
        print("\n❌ Error: query is required (unless using --check-keys)")
        sys.exit(1)

    # ── 确定输出路径 ──
    if args.output:
        output_path = _PROJECT_ROOT / "references" / args.output
    else:
        output_path = default_output_path(args.query, args.report_type)

    verbose = not args.quiet

    # ── 运行 ──
    try:
        asyncio.run(run_research(
            query=args.query,
            report_type=args.report_type,
            output_path=output_path,
            verbose=verbose,
        ))
    except KeyboardInterrupt:
        print("\n⚠️  Interrupted by user.")
        sys.exit(130)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
