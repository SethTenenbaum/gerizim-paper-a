from playwright.sync_api import sync_playwright
import sys
import time

URL = "https://whc.unesco.org/en/tentativelists/5706/"

def main():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True, args=["--no-sandbox"]) 
        context = browser.new_context(user_agent="gerizim-paper-a-archiver/1.0")
        page = context.new_page()
        print("Navigating to", URL)
        try:
            page.goto(URL, timeout=60000)
        except Exception as e:
            print("goto error:", e)
        # wait for network to settle and for JS to run
        page.wait_for_load_state("networkidle", timeout=30000)
        time.sleep(2)

        html = page.content()
        with open("unesco_5706_rendered.html", "w", encoding="utf-8") as f:
            f.write(html)
        print("Saved unesco_5706_rendered.html")

        try:
            page.screenshot(path="unesco_5706.png", full_page=True)
            print("Saved unesco_5706.png")
        except Exception as e:
            print("screenshot error:", e)

        try:
            # PDF only works in Chromium; use A4
            page.pdf(path="unesco_5706.pdf", format="A4")
            print("Saved unesco_5706.pdf")
        except Exception as e:
            print("pdf error:", e)

        browser.close()

if __name__ == "__main__":
    main()
